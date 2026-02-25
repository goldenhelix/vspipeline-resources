#!/usr/bin/env python3 

import csv
import sys
import subprocess
import json


def levenshtein_distance(s1, s2):
    """
    Calculate the Levenshtein distance between two strings.
    Returns the number of edits needed to transform s1 into s2.
    """

    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]


def get_catalog_headers():
    """
    Get the headers from the catalog.
    """

    catalog_info = json.loads(subprocess.run(["gautil", "client", "catalog-info", "SampleCatalog"], 
                                             capture_output=True, text=True).stdout)
    catalog_fields = catalog_info["schema"]["fields"]
    catalog_headers = [(f["name"], f["symbol"]) for f in catalog_fields]

    return catalog_headers


def get_manifest_headers(input_file):
    """
    Get the headers from the manifest file.
    """

    # Check file extension to determine delimiter
    if input_file.endswith('.tsv'):
        delimiter = '\t'
    else:
        delimiter = ','
        
    with open(input_file, "r") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader.fieldnames)


def get_best_match(header, headers):
    """
    Get the best match for a header from a list of headers.
    """

    best_match = None
    best_match_score = float('inf')  # Start with infinity for minimum distance
    for h in headers:
        score = levenshtein_distance(header.lower(), h[0].lower())
        if score < best_match_score:
            best_match_score = score
            best_match = h
    
    # If no good match found (threshold: distance should be reasonable relative to string length)
    if best_match is None or best_match_score > len(header) * 0.5:
        print(f"Warning: No good match found for header '{header}'", file=sys.stderr)
        return None
    
    return best_match


def match_headers_to_catalog(manifest_headers, catalog_headers):
    """
    Match the headers from the manifest file to the headers from the catalog.
    """

    matched_headers = dict()
    for header in manifest_headers:
        best_match = get_best_match(header, catalog_headers)
        matched_headers[header] = best_match
    return matched_headers


def normalize_sex_value(value):
    """
    Normalize sex values to Female, Male, or Unknown.
    """
    if not value or not value.strip():
        return "Unknown"
    
    value_lower = value.strip().lower()
    
    # Check for female variations
    if value_lower in ['f', 'female', 'woman', 'girl']:
        return "Female"
    
    # Check for male variations
    if value_lower in ['m', 'male', 'man', 'boy']:
        return "Male"
    
    # Default to Unknown for anything else
    return "Unknown"


def update_catalog(sample_line, matched_headers):
    """
    For a given sample line, update the catalog with the matched headers.
    """

    # Collect all matched headers and their values
    upsert_pairs = []
    for header, value in sample_line.items():
        if header in matched_headers and matched_headers[header] is not None:
            catalog_field = matched_headers[header][1]  # Get the field symbol from the tuple
            
            # Special handling for Sex field
            if catalog_field.lower() == 'sex':
                value = normalize_sex_value(value)
            
            upsert_pairs.append(f"{catalog_field}={value}")
    
    # Run single subprocess command with all pairs
    if upsert_pairs:
        sample_name = sample_line.get('Sample', 'Unknown Sample')
        print(f"\nUpdating catalog for sample: {sample_name}")
        print("Fields being updated:")
        for pair in upsert_pairs:
            field, value = pair.split('=', 1)
            print(f"    {field}: {value}")
        # Run the subprocess and capture output
        result = subprocess.run(
            ["gautil", "client", "catalog-upsert", "SampleCatalog"] + upsert_pairs,
            capture_output=True,
            text=True
        )
        
        # Check for errors and print output
        if result.stdout:
            print(f"Output: {result.stdout}")
        
        if result.stderr:
            print(f"Warning/Error: {result.stderr}", file=sys.stderr)
        
        if result.returncode != 0:
            print(f"Command failed with return code: {result.returncode}", file=sys.stderr)
            print(f"Failed command: gautil client catalog-upsert SampleCatalog {' '.join(upsert_pairs)}", file=sys.stderr)
        else:
            print("✓ Successfully updated catalog")

    return None


def main():
    """
    Main function to parse command line arguments and process the manifest
    """

    if len(sys.argv) < 2:
        print("Usage: python process_manifest.py <input_tsv_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    print(f"Processing manifest file: {input_file}")
    print("Getting headers from SampleCatalog...")
    catalog_headers = get_catalog_headers()
    print("Getting headers from manifest file...")
    manifest_headers = get_manifest_headers(input_file)
    print("Matching manifest headers to catalog fields...")
    matched_headers = match_headers_to_catalog(manifest_headers, catalog_headers)
    
    # Check file extension to determine delimiter
    if input_file.endswith('.tsv'):
        delimiter = '\t'
    else:
        delimiter = ','
    
    print("Updating catalog with sample data...")
    with open(input_file, "r") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for line in reader:
            update_catalog(line, matched_headers)


if __name__ == "__main__":
    main()
