#!/usr/bin/env python3

import sys
import subprocess
import requests
import json


def get_panel_genes(panel_name):
    """
    Fetch panel information from PanelApp AUS API and extract gene symbols.
    
    Args:
        panel_name: The name of the panel to fetch
        
    Returns:
        A comma-separated string of gene symbols
    """
    url = f"https://panelapp-aus.org/api/v1/panels/{panel_name}/"
    
    try:
        print(f"Fetching panel information for: {panel_name}")
        response = requests.get(url, headers={"accept": "application/json"})
        response.raise_for_status()
        
        data = response.json()
        
        # Extract gene symbols from the genes array
        genes = data.get("genes", [])
        gene_symbols = [gene["gene_data"]["gene_symbol"] for gene in genes if "gene_data" in gene and "gene_symbol" in gene["gene_data"]]
        
        gene_list = ",".join(gene_symbols)
        print(f"Found {len(gene_symbols)} genes: {gene_list}")
        
        return gene_list
        
    except requests.exceptions.RequestException as e:
        print(f"Error fetching panel data: {e}", file=sys.stderr)
        sys.exit(1)
    except (KeyError, json.JSONDecodeError) as e:
        print(f"Error parsing panel data: {e}", file=sys.stderr)
        sys.exit(1)


def update_sample_catalog(sample_name, panel_name, gene_list):
    """
    Update the SampleCatalog with Panels and GeneList fields.
    
    Args:
        sample_name: The sample name to update
        panel_name: The panel name to set
        gene_list: The comma-separated list of genes
    """
    print(f"\nUpdating catalog for sample: {sample_name}")
    print(f"  Panels: {panel_name}")
    print(f"  GeneList: {gene_list}")
    
    # Prepare the upsert command
    upsert_pairs = [
        f"Sample={sample_name}",
        f"Panels={panel_name}",
        f"GeneList={gene_list}"
    ]
    
    # Run the gautil command
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
        sys.exit(1)
    else:
        print("✓ Successfully updated catalog")


def main():
    """
    Main function to parse command line arguments and update the catalog.
    """
    if len(sys.argv) != 3:
        print("Usage: python update_catalog_with_panelapp_aus_genes.py <sample_name> <panel_name>")
        sys.exit(1)
    
    sample_name = sys.argv[1]
    panel_name = sys.argv[2]
    
    # Fetch gene list from PanelApp API
    gene_list = get_panel_genes(panel_name)
    
    # Update the catalog
    update_sample_catalog(sample_name, panel_name, gene_list)


if __name__ == "__main__":
    main()
