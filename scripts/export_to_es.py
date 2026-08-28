"""
Export Hail Table to Elasticsearch.

This script reads a Hail Table from disk and exports it to an Elasticsearch index.
It configures the Spark context required by Hail and handles the connection details.

Usage:
    uv run scripts/export_to_es.py --input-path <path> --index <index_name>
"""

import argparse
import hail as hl


def export_table(input_path: str, host: str, port: int, index: str, user: str, password: str) -> None:
    """
    Export a Hail Table to Elasticsearch.

    Args:
        input_path (str): Path to the input Hail Table directory.
        host (str): Elasticsearch host address.
        port (int): Elasticsearch port number.
        index (str): Name of the target Elasticsearch index.
    """
    # Initialize Hail with Spark configuration optimized for large datasets
    hl.init(spark_conf={
    'spark.driver.memory': '30g',
    'spark.executor.memory': '30g',
    'spark.driver.maxResultSize': '4g',
    'spark.executor.cores': '4',
    'spark.sql.shuffle.partitions': '400',
    'spark.kryoserializer.buffer.max': '512m',

    # Fix para erro de bind em ambiente Docker
    'spark.driver.bindAddress': '0.0.0.0',
    'spark.driver.host': '127.0.0.1'
    })
    
    print(f"Loading Hail Table from {input_path}...")
    ht = hl.read_table(input_path)
    
    ht = ht.naive_coalesce(3000)
    # Flatten the table structure for easier indexing in Elasticsearch
    ht = ht.flatten()
    
    # Print schema for verification
    print("Table Schema:")
    ht.describe()
    
    print(f"Exporting to Elasticsearch ({host}:{port}/{index})...")
    
    # Export to Elasticsearch
    # Note: 'es.nodes.wan.only' is set to true to allow connecting to container/remote IPs
    hl.export_elasticsearch(
        ht,
        host=host,
        port=port,
        index=index,
        index_type='_doc',
        block_size=100,
        config={                                                                                                                                         
                "es.nodes.wan.only": "true",                                                                                                                 
                "es.net.http.auth.user": user,                                                                                                               
                "es.net.http.auth.pass": password                                                                                                            
        }
    )
    
    print(f"Successfully exported table to index '{index}'.")


def get_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.

    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Export Hail Table to Elasticsearch",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input", 
        type=str, 
        required=True, 
        help="Path to the input Hail Table"
    )
    parser.add_argument(
        "--host", 
        type=str, 
        default="localhost", 
        help="Elasticsearch host"
    )
    parser.add_argument(
        "--port", 
        type=int, 
        default=9200, 
        help="Elasticsearch port"
    )
    parser.add_argument(                                                                                                                                 
            "--index",                                                                                                                                       
            type=str, 
            default="fiocruz_variants",
            help="Elasticsearch index name"
        )
    parser.add_argument(
        "--user", 
        type=str, 
        default="elastic", 
        help="Elasticsearch username"
    )
    parser.add_argument(
        "--password", 
        type=str, 
        required=True, 
        help="Elasticsearch password"
    )
    return parser


def main() -> None:
    """Main entry point."""
    parser = get_parser()
    args = parser.parse_args()
    
    export_table(
            input_path=args.input,
            host=args.host,
            port=args.port,
            index=args.index,
            user=args.user,
            password=args.password
    )


if __name__ == "__main__":
    main()
