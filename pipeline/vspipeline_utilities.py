#!/usr/bin/env python3
"""
VSPipeline Utilities

This script provides four main modes:
1. Find sample relationships - analyzes sample metadata to identify related samples
2. Generate vsbatch file - creates VSBatch files from template and sample data
3. Update CNV file - adds CNV State field to CNV files based on copy number values
4. Update vspipeline JSON - creates or updates per-sample, per-task vspipeline_inputs.json files with sample file information

Usage:
    python vspipeline_utilities.py find-relationships --samples <sample_file.txt> [--cohort]
    python vspipeline_utilities.py generate-vsbatch --template <template.vsproject-template> --samples <sample_file.txt>
    python vspipeline_utilities.py update-cnv-file --input <cnv_file.cns> --output <cnv_file_updated.cns>
    python vspipeline_utilities.py update-vspipeline-json --input_dir <directory> --sample <sample_name> --file <file_path> --file_type <multiverse|cnv|bnd|region> --task_name <task_name> --scratch_dir <scratch_directory>
"""

import argparse
import sys
import os
import csv
import json
import re
import shutil
import time
import fcntl
from pathlib import Path
from typing import List, Dict, Any, Set, Optional


def _find_axis_uuids(obj: Any, axis_lookup: Dict[str, Dict[str, Any]]) -> None:
    """
    Recursively find axisUuids in the data structure and populate the lookup dictionary.
    
    Args:
        obj: The object to search (dict, list, or primitive)
        axis_lookup: Dictionary to populate with axis UUID information
    """
    if isinstance(obj, dict):
        if 'axisUuids' in obj:
            for axis in obj['axisUuids']:
                if axis.get('className') == 'AxisUuidOptions':
                    fields = axis.get('fields', {})
                    uuid = fields.get('uuid')
                    if uuid:
                        axis_lookup[uuid] = fields
        # Recursively search nested dictionaries
        for value in obj.values():
            _find_axis_uuids(value, axis_lookup)
    elif isinstance(obj, list):
        # Recursively search list items
        for item in obj:
            _find_axis_uuids(item, axis_lookup)


def _find_table_data_options(obj: Any, path: str = "") -> List[tuple[str, Dict[str, Any]]]:
    """
    Recursively find VSTableDataOptions in the data structure.
    
    Args:
        obj: The object to search (dict, list, or primitive)
        path: Current path in the data structure (for debugging)
        
    Returns:
        List of tuples containing (path, table_option_dict)
    """
    results = []
    if isinstance(obj, dict):
        if obj.get('className') == 'VSTableDataOptions':
            results.append((path, obj))
        # Recursively search nested dictionaries
        for key, value in obj.items():
            results.extend(_find_table_data_options(value, f"{path}.{key}"))
    elif isinstance(obj, list):
        # Recursively search list items
        for i, item in enumerate(obj):
            results.extend(_find_table_data_options(item, f"{path}[{i}]"))
    return results


def _parse_template_file(template_file: str) -> Dict[str, Any]:
    """
    Parse VSProject template file to check for import algorithms and find table IDs.
    
    Args:
        template_file: Path to .vsproject-template file
        
    Returns:
        Dictionary with:
        - 'import_algorithms': Dictionary with keys 'cnv', 'breakend', 'region' and boolean values
        - 'table_ids': Dictionary with keys 'cnv', 'breakend', 'region' and table_id values (or None)
    """
    try:
        with open(template_file, 'r', encoding='utf-8') as f:
            template_data = json.load(f)
        
        # Initialize result dictionaries
        import_algorithms = {
            'cnv': False,
            'breakend': False, 
            'region': False
        }
        table_ids = {
            'cnv': None,
            'breakend': None,
            'region': None
        }
        
        # Build axis lookup dictionary
        axis_lookup = {}
        _find_axis_uuids(template_data, axis_lookup)
        
        # Find all VSTableDataOptions
        table_options = _find_table_data_options(template_data)
        
        # Collect all matching tables by axis type
        matching_tables = {
            'cnv': [],
            'breakend': [],
            'region': []
        }
        
        # Process each VSTableDataOptions to find table IDs by axis type
        for path, table_option in table_options:
            fields = table_option.get('fields', {})
            user_id = fields.get('userId', '')
            axis_uuid_list = fields.get('axisUuidList', [])
            
            # Check each axis UUID to determine the table type
            for axis_uuid in axis_uuid_list:
                if axis_uuid in axis_lookup:
                    axis_fields = axis_lookup[axis_uuid]
                    axis_type = axis_fields.get('axisType', '')
                    
                    # Collect matching tables
                    if axis_type == 'CNV':
                        matching_tables['cnv'].append(user_id)
                    elif axis_type == 'Breakend':
                        matching_tables['breakend'].append(user_id)
                    elif axis_type == 'Region':
                        matching_tables['region'].append(user_id)
        
        # Select the best table for each type using priority rules
        for table_type in ['cnv', 'breakend', 'region']:
            candidates = matching_tables[table_type]
            
            if not candidates:
                continue
            
            # Remove duplicates while preserving order
            seen = set()
            unique_candidates = []
            for candidate in candidates:
                if candidate not in seen:
                    seen.add(candidate)
                    unique_candidates.append(candidate)
            
            if len(unique_candidates) == 1:
                # Only one option
                table_ids[table_type] = unique_candidates[0]
            elif len(unique_candidates) > 1:
                # Multiple options: prefer custom names (non-"TableN") over generic ones
                custom_names = [c for c in unique_candidates if not (c.startswith('Table') and c.replace('Table', '').isdigit())]
                generic_names = [c for c in unique_candidates if c.startswith('Table') and c.replace('Table', '').isdigit()]
                
                if custom_names:
                    # Prefer custom names, alphabetically first
                    selected = sorted(custom_names)[0]
                    table_ids[table_type] = selected
                    print(f"Warning: Multiple {table_type.upper()} tables found ({', '.join(unique_candidates)}). Selected '{selected}' (preferring custom name over generic).")
                elif generic_names:
                    # Only generic names, use first alphabetically
                    selected = sorted(generic_names)[0]
                    table_ids[table_type] = selected
                    print(f"Warning: Multiple {table_type.upper()} tables found ({', '.join(unique_candidates)}). Selected '{selected}' (alphabetically first).")
                else:
                    # Fallback: use first one
                    table_ids[table_type] = unique_candidates[0]
                    print(f"Warning: Multiple {table_type.upper()} tables found ({', '.join(unique_candidates)}). Selected '{table_ids[table_type]}' (first encountered).")
        
        # Navigate to the algorithms section
        project = template_data.get('project', {})
        fields = project.get('fields', {})
        dal = fields.get('dal', {})
        dal_fields = dal.get('fields', {})
        data_trees = dal_fields.get('dataTrees', [])
        
        # Check each data tree for algorithms
        for data_tree in data_trees:
            tree_fields = data_tree.get('fields', {})
            algorithms = tree_fields.get('algorithms', [])
            
            for algorithm in algorithms:
                alg_fields = algorithm.get('fields', {})
                alg_name = alg_fields.get('algName', '')
                
                # Check for import algorithms
                if 'import' in alg_name.lower():
                    if 'cnv' in alg_name.lower():
                        import_algorithms['cnv'] = True
                    elif 'breakend' in alg_name.lower() or 'bnd' in alg_name.lower():
                        import_algorithms['breakend'] = True
                    elif 'region' in alg_name.lower():
                        import_algorithms['region'] = True
        
        return {
            'import_algorithms': import_algorithms,
            'table_ids': table_ids
        }
        
    except Exception as e:
        print(f"Error parsing template file: {e}")
        return {
            'import_algorithms': {'cnv': False, 'breakend': False, 'region': False},
            'table_ids': {'cnv': None, 'breakend': None, 'region': None}
        }


def _parse_sample_data_json(sample_data_file: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse sample data JSON file to extract file paths by type.
    
    Args:
        sample_data_file: Path to JSON file containing sample data
        
    Returns:
        Dictionary with sample names as keys and dictionaries containing file lists by type
        as values (multiverse, cnv, bnd, region)
    """
    try:
        with open(sample_data_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        samples = data.get('samples', {})
        result = {}
        
        for sample_name, sample_info in samples.items():
            result[sample_name] = {
                'multiverse': sample_info.get('multiverse', []),
                'cnv': sample_info.get('cnv', []),
                'bnd': sample_info.get('bnd', []),
                'region': sample_info.get('region', [])
            }
        
        return result
        
    except Exception as e:
        print(f"Error parsing sample data JSON: {e}")
        return {}


def _merge_vspipeline_json_files(json_files: List[str]) -> Dict[str, Dict[str, List[str]]]:
    """
    Merge multiple per-sample, per-task vspipeline_inputs.json files into a single structure.
    
    Args:
        json_files: List of paths to JSON files to merge
        
    Returns:
        Dictionary with sample names as keys and dictionaries containing file lists by type
        as values (multiverse, cnv, bnd, region)
    """
    merged_data = {}
    
    for json_file in json_files:
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            samples = data.get('samples', {})
            
            for sample_name, sample_info in samples.items():
                if sample_name not in merged_data:
                    merged_data[sample_name] = {
                        'multiverse': [],
                        'cnv': [],
                        'bnd': [],
                        'region': []
                    }
                
                # Merge file lists, avoiding duplicates
                for file_type in ['multiverse', 'cnv', 'bnd', 'region']:
                    files = sample_info.get(file_type, [])
                    for file_path in files:
                        if file_path not in merged_data[sample_name][file_type]:
                            merged_data[sample_name][file_type].append(file_path)
        
        except Exception as e:
            print(f"Warning: Error parsing JSON file {json_file}: {e}")
            continue
    
    return merged_data


def _generate_vsbatch_content(template_file: str, 
                              sample_group: List[str], 
                              sample_data: Dict[str, Dict[str, List[str]]], 
                              import_algorithms: Dict[str, bool],
                              table_ids: Dict[str, Optional[str]],
                              overwrite_project: bool = True) -> str:
    """
    Generate VSBatch file content for a single group of related samples.
    
    Args:
        template_file: Path to the template file
        sample_group: List of sample names in this group
        sample_data: Parsed sample data with file paths by type
        import_algorithms: Dictionary indicating which import algorithms are available
        table_ids: Dictionary with table_id values for 'cnv', 'breakend', and 'region' (may be None)
        overwrite_project: Whether to overwrite existing projects (default: True)
        
    Returns:
        VSBatch file content as a string
    """
    # Determine project name based on number of samples
    if len(sample_group) < 6:
        # For small groups, use underscore-separated list of samples
        project_name = '_'.join(sorted(sample_group))
    else:
        # For large groups, use first sample name followed by _cohort
        project_name = f"{sample_group[0]}_cohort"
    
    vsbatch_lines = []
    
    # Add project creation command
    template_name = os.path.basename(template_file)
    project_create_cmd = f'project_create "AppData/Projects/{project_name}" template="{template_name}"'
    if overwrite_project:
        project_create_cmd += ' overwrite=true'
    vsbatch_lines.append(project_create_cmd)
    vsbatch_lines.append('download_required_sources')
    
    # Collect all files for this sample group
    multiverse_files = []
    cnv_files = []
    bnd_files = []
    region_files = []
    
    for sample_name in sample_group:
        if sample_name in sample_data:
            sample_info = sample_data[sample_name]
            multiverse_files.extend(sample_info.get('multiverse', []))
            cnv_files.extend(sample_info.get('cnv', []))
            bnd_files.extend(sample_info.get('bnd', []))
            region_files.extend(sample_info.get('region', []))
    
    # Add multiverse import (always present)
    if multiverse_files:
        # Quote each file and join with commas (no spaces)
        quoted_files = ','.join([f'"{f}"' for f in multiverse_files])
        vsbatch_lines.append(f'import files={quoted_files} \
            sample_fields_catalog=SampleCatalog \
            append_files_with_same_samples=True')
        vsbatch_lines.append(f'task_wait')
    
    # Add CNV import if algorithm exists and files are available
    if import_algorithms.get('cnv', False) and cnv_files:
        quoted_files = ','.join([f'"{f}"' for f in cnv_files])
        table_id = table_ids.get('cnv') or 'Table1'  # Fallback to 'Table1' if None or not found
        vsbatch_lines.append(f'update_cnv_import files={quoted_files} table_id="{table_id}"')
        vsbatch_lines.append(f'task_wait')

    # Add breakend import if algorithm exists and files are available
    if import_algorithms.get('breakend', False) and bnd_files:
        quoted_files = ','.join([f'"{f}"' for f in bnd_files])
        table_id = table_ids.get('breakend') or 'Table1'  # Fallback to 'Table1' if None or not found
        vsbatch_lines.append(f'update_bnd_import files={quoted_files} table_id="{table_id}"')
        vsbatch_lines.append(f'task_wait')
    
    # Add region import if algorithm exists and files are available
    if import_algorithms.get('region', False) and region_files:
        quoted_files = ','.join([f'"{f}"' for f in region_files])
        table_id = table_ids.get('region') or 'Table1'  # Fallback to 'Table1' if None or not found
        vsbatch_lines.append(f'update_region_import files={quoted_files} table_id="{table_id}"')
        vsbatch_lines.append(f'task_wait')
    
    # Add standard export and project management commands
    vsbatch_lines.extend([
        'get_task_list',
        'project_save',
        'project_close'
    ])
    
    return '\n'.join(vsbatch_lines)


def _parse_sample_file(sample_file: str) -> tuple[Dict[str, List[str]], Set[str]]:
    """Parse TSV file and return sample data and all sample names."""
    sample_data = {}
    all_samples = set()
    
    with open(sample_file, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip header row
        
        for row in reader:
            if not row:
                continue
            
            sample_name = row[0].strip()
            if not sample_name:
                continue
            
            all_samples.add(sample_name)
            sample_data[sample_name] = [value.strip() for value in row[1:] if value and value.strip()]
    
    return sample_data, all_samples


def _find_connected_samples(sample: str, sample_data: Dict[str, List[str]], 
                          all_samples: Set[str], processed_samples: Set[str]) -> Set[str]:
    """Find all samples connected to the given sample through relationships."""
    current_group = {sample}
    processed_samples.add(sample)
    to_process = [sample]
    
    while to_process:
        current_sample = to_process.pop(0)
        
        # Check if current_sample appears in any other sample's metadata
        for other_sample, metadata_values in sample_data.items():
            if other_sample in processed_samples:
                continue
            
            if current_sample in metadata_values:
                current_group.add(other_sample)
                processed_samples.add(other_sample)
                to_process.append(other_sample)
        
        # Check if any sample in current_group appears in current_sample's metadata
        for metadata_value in sample_data.get(current_sample, []):
            if metadata_value in all_samples and metadata_value not in processed_samples:
                current_group.add(metadata_value)
                processed_samples.add(metadata_value)
                to_process.append(metadata_value)
    
    return current_group


def _output_sample_groups(relationships: Dict[str, List[str]], output_file: str) -> None:
    """Output sample groups to CSV file."""
    print(f"Found {len(relationships)} sample groups:")
    
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        for group_key in sorted(relationships.keys()):
            group = relationships[group_key]
            print(f"Group: {', '.join(group)}")
            f.write(','.join(group) + '\n')
    
    print(f"\nTotal groups found: {len(relationships)}")
    print(f"Results written to: {output_file}")


def find_sample_relationships(sample_file: str, output_file: str, cohort: bool = False) -> None:
    """
    Find sample relationships based on metadata in the provided sample file.
    
    Args:
        sample_file: Path to a TSV file containing sample metadata
        output_file: Path to output CSV file for sample groups
        cohort: If True, put all samples into a single cohort group instead of finding relationships
    """
    if cohort:
        print(f"Creating cohort group from: {sample_file}")
    else:
        print(f"Finding sample relationships from: {sample_file}")
    
    if not os.path.exists(sample_file):
        print(f"Error: Sample file '{sample_file}' does not exist.")
        sys.exit(1)
    
    try:
        sample_data, all_samples = _parse_sample_file(sample_file)
        print(f"Found {len(all_samples)} samples")
        
        relationships = {}
        
        if cohort:
            # Put all samples into a single cohort group
            if all_samples:
                relationships['cohort'] = sorted(all_samples)
                print(f"Created single cohort group with {len(all_samples)} samples")
        else:
            # Find relationships between samples
            processed_samples = set()
            
            for sample in all_samples:
                if sample in processed_samples:
                    continue
                
                current_group = _find_connected_samples(sample, sample_data, all_samples, processed_samples)
                group_key = min(current_group)
                relationships[group_key] = sorted(current_group)
        
        _output_sample_groups(relationships, output_file)
        
    except Exception as e:
        print(f"Error processing sample file: {e}")
        sys.exit(1)


def update_cnv_file(input_file: str, output_file: str) -> None:
    """
    Update CNV file by adding a CNV State field based on copy number values.
    Supports both 'cn' (integer) and 'copynumber_state' (decimal) fields.
    
    Args:
        input_file: Path to input CNV file (TSV format)
        output_file: Path to output CNV file with CNV State field added
    """
    print(f"Updating CNV file: {input_file}")
    print(f"Output file: {output_file}")
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        sys.exit(1)
    
    try:
        # Read the input file
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        if not lines:
            print("Error: Input file is empty.")
            sys.exit(1)
        
        # Find the header line (look for tab-separated column names)
        header_line = None
        header_line_num = 0
        
        # First, try to find a line that looks like column headers
        for i, line in enumerate(lines):
            line = line.strip()
            if line and '\t' in line:
                # Check if this line contains typical column names
                headers_candidate = line.split('\t')
                if any(col in ['chr', 'start', 'end', 'cn', 'copynumber_state', 'coverage', 'confidence'] for col in headers_candidate):
                    header_line = line
                    header_line_num = i
                    break
        
        if not header_line:
            print("Error: No header line found in input file.")
            sys.exit(1)
        
        headers = header_line.split('\t')
        
        # Check for copy number column (cn or copynumber_state)
        cn_index = None
        cn_field_name = None
        
        if 'cn' in headers:
            cn_index = headers.index('cn')
            cn_field_name = 'cn'
        elif 'copynumber_state' in headers:
            cn_index = headers.index('copynumber_state')
            cn_field_name = 'copynumber_state'
        else:
            print("Error: Neither 'cn' nor 'copynumber_state' column found in input file.")
            print(f"Available columns: {', '.join(headers)}")
            sys.exit(1)
        
        print(f"Using copy number field: '{cn_field_name}'")
        
        # Add CNV State column to header
        headers.append('CNV State')
        updated_lines = ['\t'.join(headers) + '\n']
        
        # Process data lines (skip comment lines and header)
        processed_count = 0
        for line_num, line in enumerate(lines[header_line_num + 1:], header_line_num + 2):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) <= cn_index:
                print(f"Warning: Line {line_num} has insufficient columns, skipping.")
                continue
            
            try:
                # Handle both integer and decimal copy number values
                cn_value = float(fields[cn_index])
                
                # Apply CNV State logic (round to nearest integer for comparison)
                cn_int = round(cn_value)
                
                if cn_int == 2:
                    cnv_state = "Normal"
                elif cn_int == 0:
                    cnv_state = "Deletion"
                elif cn_int == 1:
                    cnv_state = "Het Deletion"
                elif cn_int >= 3:
                    cnv_state = "Duplicate"
                else:
                    cnv_state = "?"
                
                # Add CNV State to the end of the line
                fields.append(cnv_state)
                updated_lines.append('\t'.join(fields) + '\n')
                processed_count += 1
                
            except ValueError:
                print(f"Warning: Line {line_num} has non-numeric copy number value '{fields[cn_index]}', skipping.")
                continue
        
        # Write the updated file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.writelines(updated_lines)
        
        print(f"Successfully processed {processed_count} CNV entries")
        print(f"Updated file written to: {output_file}")
        
        # Show a preview of the first few lines
        print("\nPreview of updated file:")
        print("=" * 50)
        with open(output_file, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i < 5:  # Show first 5 lines
                    print(line.strip())
                else:
                    break
        print("=" * 50)
        
    except Exception as e:
        print(f"Error updating CNV file: {e}")
        sys.exit(1)


def _safe_read_json_file(file_path: str, max_retries: int = 30, retry_delay: float = 0.5) -> Dict[str, Any]:
    """
    Safely read a JSON file with retry logic for file locking.
    
    Args:
        file_path: Path to the JSON file
        max_retries: Maximum number of retry attempts
        retry_delay: Initial delay between retries (increases exponentially)
        
    Returns:
        Dictionary containing the JSON data, or empty dict if file doesn't exist
    """
    for attempt in range(max_retries):
        try:
            if not os.path.exists(file_path):
                return {}
            
            # Try to open and lock the file
            with open(file_path, 'r', encoding='utf-8') as f:
                try:
                    # Try to acquire a shared lock (non-blocking) for reading
                    fcntl.flock(f.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
                    
                    # Read the file content
                    content = f.read()
                    
                    if not content.strip():
                        return {}
                    
                    return json.loads(content)
                    
                except (IOError, OSError):
                    # File is locked by another process
                    if attempt < max_retries - 1:
                        # Exponential backoff
                        delay = retry_delay * (2 ** min(attempt, 5))
                        time.sleep(delay)
                        continue
                    else:
                        raise Exception(f"Could not acquire lock on {file_path} after {max_retries} attempts")
                        
        except json.JSONDecodeError as e:
            # If JSON is invalid, wait a bit and retry (might be mid-write)
            if attempt < max_retries - 1:
                delay = retry_delay * (2 ** min(attempt, 5))
                time.sleep(delay)
                continue
            else:
                raise Exception(f"Invalid JSON in {file_path}: {e}")
        except Exception as e:
            if attempt < max_retries - 1:
                delay = retry_delay * (2 ** min(attempt, 5))
                time.sleep(delay)
                continue
            else:
                raise
    
    return {}


def _safe_write_json_file(file_path: str, data: Dict[str, Any], max_retries: int = 30, retry_delay: float = 0.5) -> None:
    """
    Safely write a JSON file with retry logic for file locking.
    
    Args:
        file_path: Path to the JSON file
        data: Dictionary to write as JSON
        max_retries: Maximum number of retry attempts
        retry_delay: Initial delay between retries (increases exponentially)
    """
    # Ensure directory exists (if there's a directory component)
    dir_path = os.path.dirname(file_path)
    if dir_path:
        os.makedirs(dir_path, exist_ok=True)
    
    for attempt in range(max_retries):
        try:
            # Open file in read-write mode
            with open(file_path, 'r+', encoding='utf-8') as f:
                try:
                    # Try to acquire an exclusive lock (non-blocking)
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    
                    # Write the JSON data
                    f.seek(0)
                    json.dump(data, f, indent=2, ensure_ascii=False)
                    f.truncate()
                    f.flush()
                    os.fsync(f.fileno())
                    
                    return  # Success
                    
                except (IOError, OSError):
                    # File is locked by another process
                    if attempt < max_retries - 1:
                        # Exponential backoff
                        delay = retry_delay * (2 ** min(attempt, 5))
                        time.sleep(delay)
                        continue
                    else:
                        raise Exception(f"Could not acquire lock on {file_path} after {max_retries} attempts")
                        
        except FileNotFoundError:
            # File doesn't exist, create it
            with open(file_path, 'w', encoding='utf-8') as f:
                try:
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    json.dump(data, f, indent=2, ensure_ascii=False)
                    f.flush()
                    os.fsync(f.fileno())
                    return
                except (IOError, OSError):
                    if attempt < max_retries - 1:
                        delay = retry_delay * (2 ** min(attempt, 5))
                        time.sleep(delay)
                        continue
                    else:
                        raise Exception(f"Could not acquire lock on {file_path} after {max_retries} attempts")
        except Exception as e:
            if attempt < max_retries - 1:
                delay = retry_delay * (2 ** min(attempt, 5))
                time.sleep(delay)
                continue
            else:
                raise


def update_vspipeline_json(input_dir: str, sample_name: str, file_path: str, file_type: str, task_name: str, scratch_dir: str) -> None:
    """
    Create or update per-sample, per-task vspipeline_inputs.json file with sample file information.
    Uses a scratch directory to avoid race conditions by copying the file to scratch, editing it,
    and then copying it back.
    
    Args:
        input_dir: Directory containing the vspipeline_inputs.json file (output directory)
        sample_name: Name of the sample
        file_path: Path to the file to add
        file_type: Type of file ('multiverse', 'cnv', 'bnd', or 'region')
        task_name: Name of the task (e.g., 'hiphase', 'methbat', 'mitorsaw')
        scratch_dir: Scratch directory for temporary file operations (e.g., '/scratch')
    """
    # Validate file type
    valid_file_types = ['multiverse', 'cnv', 'bnd', 'region']
    if file_type not in valid_file_types:
        print(f"Error: Invalid file type '{file_type}'. Must be one of: {', '.join(valid_file_types)}")
        sys.exit(1)
    
    # Construct file names
    json_filename = f'{sample_name}_{task_name}_vspipeline_inputs.json'
    output_json_path = os.path.join(input_dir, json_filename)
    scratch_json_path = os.path.join(scratch_dir, json_filename)
    
    print(f"Updating vspipeline_inputs.json in: {input_dir}")
    print(f"Sample: {sample_name}")
    print(f"Task: {task_name}")
    print(f"File: {file_path}")
    print(f"Type: {file_type}")
    print(f"Using scratch directory: {scratch_dir}")
    
    try:
        # Ensure scratch directory exists
        os.makedirs(scratch_dir, exist_ok=True)
        
        # Copy existing JSON file from output directory to scratch (if it exists)
        if os.path.exists(output_json_path):
            print(f"Copying existing JSON file from {output_json_path} to {scratch_json_path}")
            shutil.copy2(output_json_path, scratch_json_path)
        else:
            print(f"JSON file does not exist yet, will create new one in scratch")
        
        # Read existing JSON file from scratch (or get empty dict if it doesn't exist)
        data = _safe_read_json_file(scratch_json_path)
        
        # Initialize structure if needed
        if 'samples' not in data:
            data['samples'] = {}
        
        # Initialize sample entry if needed
        if sample_name not in data['samples']:
            data['samples'][sample_name] = {
                'multiverse': [],
                'cnv': [],
                'bnd': [],
                'region': []
            }
        
        # Add file to the appropriate list (avoid duplicates)
        file_list = data['samples'][sample_name][file_type]
        if file_path not in file_list:
            file_list.append(file_path)
            print(f"Added file to {file_type} list for sample {sample_name}")
        else:
            print(f"File already exists in {file_type} list for sample {sample_name}")
        
        # Write to scratch directory
        _safe_write_json_file(scratch_json_path, data)
        print(f"Successfully updated JSON in scratch: {scratch_json_path}")
        
        # Copy updated file back to output directory
        print(f"Copying updated JSON file from {scratch_json_path} to {output_json_path}")
        shutil.copy2(scratch_json_path, output_json_path)
        
        print(f"Successfully updated: {output_json_path}")
        
        # Show current state for this sample
        print(f"\nCurrent state for sample '{sample_name}':")
        sample_data = data['samples'][sample_name]
        for ft in valid_file_types:
            count = len(sample_data.get(ft, []))
            print(f"  {ft}: {count} file(s)")
        
    except Exception as e:
        print(f"Error updating vspipeline_inputs.json: {e}")
        sys.exit(1)


def generate_vsbatch_file(template_file: str, sample_file: str, sample_data_files: List[str], overwrite_project: bool = True, scratch_dir: str = None) -> None:
    """
    Generate VSBatch file from template and sample data.
    Uses a scratch directory to copy JSON files before reading them.
    
    Args:
        template_file: Path to .vsproject-template file
        sample_file: Path to CSV file containing sample groups
        sample_data_files: List of paths to JSON files containing sample data with file paths (will be merged)
        overwrite_project: Whether to overwrite existing projects (default: True)
        scratch_dir: Optional scratch directory for copying JSON files before reading (e.g., '/scratch')
    """
    print(f"Generating VSBatch file from template: {template_file}")
    print(f"Using sample groups from: {sample_file}")
    print(f"Using sample data from {len(sample_data_files)} JSON file(s)")
    if scratch_dir:
        print(f"Using scratch directory: {scratch_dir}")
    
    # Validate files exist
    if not os.path.exists(template_file):
        print(f"Error: Template file '{template_file}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(sample_file):
        print(f"Error: Sample file '{sample_file}' does not exist.")
        sys.exit(1)
    
    # Validate all sample data files exist
    for sample_data_file in sample_data_files:
        if not os.path.exists(sample_data_file):
            print(f"Warning: Sample data file '{sample_data_file}' does not exist, skipping.")
    
    # Filter to only existing files
    existing_files = [f for f in sample_data_files if os.path.exists(f)]
    if not existing_files:
        print(f"Error: No valid sample data files found.")
        sys.exit(1)
    
    try:
        # If scratch directory is provided, copy JSON files there before reading
        files_to_read = existing_files
        if scratch_dir:
            # Ensure scratch directory exists
            os.makedirs(scratch_dir, exist_ok=True)
            print(f"Copying {len(existing_files)} JSON file(s) to scratch directory...")
            scratch_files = []
            for json_file in existing_files:
                json_filename = os.path.basename(json_file)
                scratch_file = os.path.join(scratch_dir, json_filename)
                shutil.copy2(json_file, scratch_file)
                scratch_files.append(scratch_file)
                print(f"  Copied {json_filename} to {scratch_file}")
            files_to_read = scratch_files
        
        # Parse template to check for import algorithms and find table IDs
        print("Parsing template file...")
        template_info = _parse_template_file(template_file)
        import_algorithms = template_info['import_algorithms']
        table_ids = template_info['table_ids']
        print(f"Available import algorithms: {import_algorithms}")
        print(f"Found table IDs: {table_ids}")
        
        # Merge sample data from all JSON files
        print(f"Merging sample data from {len(files_to_read)} JSON file(s)...")
        sample_data = _merge_vspipeline_json_files(files_to_read)
        print(f"Found data for {len(sample_data)} samples after merging")
        
        # Parse sample groups CSV
        print("Parsing sample groups...")
        sample_groups = []
        with open(sample_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:
                    sample_groups.append([s.strip() for s in line.split(',')])
        
        print(f"Found {len(sample_groups)} sample groups")
        
        # Generate separate VSBatch files for each group
        print("Generating VSBatch files...")
        generated_files = []
        
        for group_idx, sample_group in enumerate(sample_groups, 1):
            group_name = f"Group_{group_idx}"
            print(f"Processing sample group {group_idx}: {', '.join(sample_group)}")
            
            # Generate VSBatch content for this group
            vsbatch_content = _generate_vsbatch_content(template_file, sample_group, sample_data, import_algorithms, table_ids, overwrite_project)
            
            # Write VSBatch file
            output_file = f"batch{group_idx}.vsbatch"
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(vsbatch_content)
            
            generated_files.append(output_file)
            print(f"Generated: {output_file}")
        
        print(f"\nGenerated {len(generated_files)} VSBatch files:")
        for file in generated_files:
            print(f"  • {file}")
        
        # Show preview of first file
        if generated_files:
            print(f"\nContent preview of {generated_files[0]}:")
            print("=" * 50)
            with open(generated_files[0], 'r', encoding='utf-8') as f:
                print(f.read())
            print("=" * 50)
        
    except Exception as e:
        print(f"Error generating VSBatch file: {e}")
        sys.exit(1)


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="VSPipeline utilities for PacBio WDL Somatic workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Find sample relationships
  python vspipeline_utilities.py find-relationships --samples samples.tsv --output sample_groups.csv
  
  # Create a single cohort group from all samples
  python vspipeline_utilities.py find-relationships --samples samples.tsv --output sample_groups.csv --cohort
  
  # Generate VSBatch files
  python vspipeline_utilities.py generate-vsbatch --template template.vsproject-template --samples sample_groups.csv --sample_files sample1_hiphase_vspipeline_inputs.json sample2_methbat_vspipeline_inputs.json
  
  # Update CNV file with CNV State field
  python vspipeline_utilities.py update-cnv-file --input cnv_file.cns --output cnv_file_updated.cns
  
  # Update vspipeline_inputs.json file
  python vspipeline_utilities.py update-vspipeline-json --input_dir /path/to/dir --sample Sample1 --file /path/to/file.vcf.gz --file_type multiverse --task_name hiphase --scratch_dir /scratch
        """
    )
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(dest='mode', help='Available modes')
    
    # Find relationships mode
    find_parser = subparsers.add_parser(
        'find-relationships',
        help='Find sample relationships from metadata file'
    )
    find_parser.add_argument(
        '--samples',
        required=True,
        help='Path to TSV file containing sample metadata'
    )
    find_parser.add_argument(
        '--output',
        required=True,
        help='Path to output CSV file for sample groups'
    )
    find_parser.add_argument(
        '--cohort',
        action='store_true',
        help='Put all samples into a single cohort group instead of finding relationships'
    )
    
    # Generate VSBatch mode
    generate_parser = subparsers.add_parser(
        'generate-vsbatch',
        help='Generate VSBatch file from template and sample data'
    )
    generate_parser.add_argument(
        '--template',
        required=True,
        help='Path to .vsproject-template file'
    )
    generate_parser.add_argument(
        '--samples',
        required=True,
        help='Path to CSV file containing sample groups'
    )
    generate_parser.add_argument(
        '--sample_files',
        required=True,
        nargs='+',
        help='Path(s) to JSON file(s) containing sample data with file paths (multiple files will be merged)'
    )
    generate_parser.add_argument(
        '--overwrite_project',
        type=str,
        default='true',
        help='Whether to overwrite existing projects (true/false, default: true)'
    )
    generate_parser.add_argument(
        '--scratch_dir',
        type=str,
        default=None,
        help='Optional scratch directory for copying JSON files before reading (e.g., /scratch)'
    )
    
    # Update CNV file mode
    update_cnv_parser = subparsers.add_parser(
        'update-cnv-file',
        help='Update CNV file by adding CNV State field based on copy number values'
    )
    update_cnv_parser.add_argument(
        '--input',
        required=True,
        help='Path to input CNV file (TSV format)'
    )
    update_cnv_parser.add_argument(
        '--output',
        required=True,
        help='Path to output CNV file with CNV State field added'
    )
    
    # Update vspipeline JSON mode
    update_json_parser = subparsers.add_parser(
        'update-vspipeline-json',
        help='Create or update per-sample, per-task vspipeline_inputs.json file with sample file information'
    )
    update_json_parser.add_argument(
        '--input_dir',
        required=True,
        help='Directory containing the vspipeline_inputs.json file'
    )
    update_json_parser.add_argument(
        '--sample',
        required=True,
        help='Name of the sample'
    )
    update_json_parser.add_argument(
        '--file',
        required=True,
        help='Path to the file to add'
    )
    update_json_parser.add_argument(
        '--file_type',
        required=True,
        choices=['multiverse', 'cnv', 'bnd', 'region'],
        help='Type of file: multiverse, cnv, bnd, or region'
    )
    update_json_parser.add_argument(
        '--task_name',
        required=True,
        help='Name of the task (e.g., hiphase, methbat, mitorsaw)'
    )
    update_json_parser.add_argument(
        '--scratch_dir',
        required=True,
        help='Scratch directory for temporary file operations (e.g., /scratch)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if no mode was specified
    if not args.mode:
        parser.print_help()
        sys.exit(1)
    
    # Execute the appropriate mode
    if args.mode == 'find-relationships':
        find_sample_relationships(args.samples, args.output, cohort=args.cohort)
    elif args.mode == 'generate-vsbatch':
        # Convert string to boolean
        overwrite = args.overwrite_project.lower() in ('true', '1', 'yes', 'on')
        # args.sample_files is already a list due to nargs='+'
        generate_vsbatch_file(args.template, args.samples, args.sample_files, overwrite_project=overwrite, scratch_dir=args.scratch_dir)
    elif args.mode == 'update-cnv-file':
        update_cnv_file(args.input, args.output)
    elif args.mode == 'update-vspipeline-json':
        update_vspipeline_json(args.input_dir, args.sample, args.file, args.file_type, args.task_name, args.scratch_dir)
    else:
        print(f"Error: Unknown mode '{args.mode}'")
        sys.exit(1)


if __name__ == '__main__':
    main()
