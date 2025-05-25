"""
VCF Analysis Agent CLI Entrypoint

- Provides a command-line interface for interacting with the agent
- Accepts prompt strings for tool invocation (e.g., validation, echo)
- Supports mock response for testing
- Supports multiple model providers (Ollama, OpenAI, Cerebras)
- Provides LanceDB integration commands:
    * init-lancedb: Initialize LanceDB table
    * add-variant: Add a variant record
    * search-embedding: Search by embedding vector

Usage examples:
    $ python -m vcf_agent.cli init-lancedb
    $ python -m vcf_agent.cli add-variant --variant_id rs123 --chrom 1 --pos 12345 --ref A --alt G --embedding 0.1,0.2,...
    $ python -m vcf_agent.cli search-embedding --embedding 0.1,0.2,...
"""

import argparse
import os
from typing import Literal, cast
import json
import sys

def main():
    """
    Command-line interface for the VCF Analysis Agent.

    Accepts a prompt string to invoke agent tools (e.g., validation, echo).
    Supports mock response for testing via the VCF_AGENT_CLI_MOCK_RESPONSE environment variable.
    Supports --raw / --no-think flag to disable chain-of-thought reasoning (raw output mode).
    Supports --model to select the model provider (ollama, openai, cerebras).
    Provides LanceDB integration commands.
    """
    mock_response = os.environ.get("VCF_AGENT_CLI_MOCK_RESPONSE")
    if mock_response is not None:
        print(mock_response)
        return
        
    parser = argparse.ArgumentParser(description="VCF Analysis Agent CLI")
    subparsers = parser.add_subparsers(dest="command")

    # Default: agent prompt
    parser.add_argument("prompt", type=str, nargs="?", help="Prompt for the agent")
    parser.add_argument(
        "--raw", "--no-think", action="store_true", 
        help="Disable chain-of-thought reasoning (raw output mode)"
    )
    parser.add_argument(
        "--model", type=str, choices=["ollama", "openai", "cerebras"], default="ollama",
        help="Select model provider (default: ollama)"
    )
    parser.add_argument(
        "--credentials", type=str, 
        help="Path to JSON credentials file for API access"
    )
    parser.add_argument(
        "--reference", "--fasta", type=str, default=None,
        help="Path to reference FASTA for normalization (required for some tools)"
    )

    # LanceDB: init
    parser_init = subparsers.add_parser("init-lancedb", help="Initialize LanceDB table for variants")
    parser_init.add_argument("--db_path", type=str, default="./lancedb", help="Path to LanceDB directory")
    parser_init.add_argument("--table_name", type=str, default="variants", help="Table name")

    # LanceDB: add variant
    parser_add = subparsers.add_parser("add-variant", help="Add a variant record to LanceDB")
    parser_add.add_argument("--db_path", type=str, default="./lancedb")
    parser_add.add_argument("--table_name", type=str, default="variants")
    parser_add.add_argument("--variant_id", type=str, required=True)
    parser_add.add_argument("--chrom", type=str, required=True)
    parser_add.add_argument("--pos", type=int, required=True)
    parser_add.add_argument("--ref", type=str, required=True)
    parser_add.add_argument("--alt", type=str, required=True)
    parser_add.add_argument("--embedding", type=str, required=True, help="Comma-separated embedding vector")
    parser_add.add_argument("--clinical_significance", type=str, default=None)

    # LanceDB: search by embedding
    parser_search = subparsers.add_parser("search-embedding", help="Search variants by embedding vector")
    parser_search.add_argument("--db_path", type=str, default="./lancedb")
    parser_search.add_argument("--table_name", type=str, default="variants")
    parser_search.add_argument("--embedding", type=str, required=True, help="Comma-separated embedding vector")
    parser_search.add_argument("--limit", type=int, default=10)
    parser_search.add_argument("--filter_sql", type=str, default=None, help="Optional SQL filter string")

    # LanceDB: update variant
    parser_update = subparsers.add_parser("update-variant", help="Update a variant record in LanceDB")
    parser_update.add_argument("--db_path", type=str, default="./lancedb")
    parser_update.add_argument("--table_name", type=str, default="variants")
    parser_update.add_argument("--variant_id", type=str, required=True, help="ID of the variant to update")
    parser_update.add_argument("--updates", type=str, required=True, help="JSON string of updates, e.g., '{\"clinical_significance\": \"Benign\", \"pos\": 12346}'")

    # LanceDB: delete variants
    parser_delete = subparsers.add_parser("delete-variants", help="Delete variants from LanceDB based on a filter")
    parser_delete.add_argument("--db_path", type=str, default="./lancedb")
    parser_delete.add_argument("--table_name", type=str, default="variants")
    parser_delete.add_argument("--filter_sql", type=str, required=True, help="SQL filter string for deletion")

    # LanceDB: create scalar index
    parser_index = subparsers.add_parser("create-lancedb-index", help="Create a scalar index on a table column")
    parser_index.add_argument("--db_path", type=str, default="./lancedb")
    parser_index.add_argument("--table_name", type=str, default="variants")
    parser_index.add_argument("--column", type=str, required=True, help="Name of the column to index")
    parser_index.add_argument("--index_type", type=str, default=None, help="Optional: Type of index (e.g., BTREE, BITMAP)")
    parser_index.add_argument("--replace", action="store_true", help="Replace existing index if it exists")

    # LanceDB: filter variants by metadata
    parser_filter = subparsers.add_parser("filter-lancedb", help="Filter variants based on metadata using a SQL query.")
    parser_filter.add_argument("--db_path", type=str, default="./lancedb")
    parser_filter.add_argument("--table_name", type=str, default="variants")
    parser_filter.add_argument("--filter_sql", type=str, required=True, help="SQL filter string for metadata query (e.g., \"chrom = '1' AND pos > 1000\")")
    parser_filter.add_argument("--select_columns", type=str, default=None, help="Comma-separated list of columns to return (e.g., variant_id,chrom,pos)")
    parser_filter.add_argument("--limit", type=int, default=None, help="Maximum number of records to return")

    args = parser.parse_args()

    # LanceDB integration commands
    if args.command == "init-lancedb":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        print(f"Initialized LanceDB table '{args.table_name}' at '{args.db_path}'")
        return
    elif args.command == "add-variant":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table, add_variants
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        embedding = [float(x) for x in args.embedding.split(",")]
        variant = {
            "variant_id": args.variant_id,
            "chrom": args.chrom,
            "pos": args.pos,
            "ref": args.ref,
            "alt": args.alt,
            "embedding": embedding,
            "clinical_significance": args.clinical_significance,
        }
        add_variants(table, [variant])
        print(f"Added variant {args.variant_id} to LanceDB.")
        return
    elif args.command == "search-embedding":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table, search_by_embedding
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        embedding = [float(x) for x in args.embedding.split(",")]
        results = search_by_embedding(table, embedding, limit=args.limit, filter_sql=args.filter_sql)
        print(results)
        return
    elif args.command == "update-variant":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table, update_variant
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        try:
            updates = json.loads(args.updates)
            if not isinstance(updates, dict):
                raise ValueError("Updates must be a valid JSON object (dictionary).")
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON string for updates: {e}", file=sys.stderr)
            sys.exit(1)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        
        try:
            update_variant(table, args.variant_id, updates)
            print(f"Update command processed for variant {args.variant_id}.")
        except Exception as e:
            print(f"Error during variant update: {e}", file=sys.stderr)
            sys.exit(1)
        return
    elif args.command == "delete-variants":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table, delete_variants
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        delete_variants(table, args.filter_sql)
        print(f"Delete command processed with filter: {args.filter_sql}")
        return
    elif args.command == "create-lancedb-index":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table, create_scalar_index
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        # Basic kwargs handling for now, could be expanded if specific index params are common
        # For example, if num_partitions was a CLI arg, it would be passed here.
        index_kwargs = {}
        create_scalar_index(table, args.column, index_type=args.index_type, replace=args.replace, **index_kwargs)
        print(f"Index creation command processed for column '{args.column}' on table '{args.table_name}'.")
        return
    elif args.command == "filter-lancedb":
        from vcf_agent.lancedb_integration import get_db, get_or_create_table, filter_variants_by_metadata
        db = get_db(args.db_path)
        table = get_or_create_table(db, args.table_name)
        select_cols_list = args.select_columns.split(',') if args.select_columns else None
        results_df = filter_variants_by_metadata(table, args.filter_sql, select_columns=select_cols_list, limit=args.limit)
        print(f"Found {len(results_df)} variants matching filter:")
        print(results_df.to_string())
        return

    # Default: agent prompt
    if args.raw:
        os.environ["VCF_AGENT_RAW_MODE"] = "1"

    # Import after setting env vars
    from vcf_agent.agent import get_agent_with_session
    from vcf_agent.config import SessionConfig
    
    # Create session config with CLI options
    model_provider = cast(Literal["ollama", "openai", "cerebras"], args.model)
    session_config = SessionConfig(
        raw_mode=args.raw if args.raw else None,
        model_provider=model_provider,
        credentials_file=args.credentials,
        reference_fasta=args.reference
    )
    
    # Get agent with specified model
    agent = get_agent_with_session(
        session_config=session_config,
        model_provider=model_provider
    )
    
    # Run prompt
    if args.prompt:
        response = agent(args.prompt)
        print(response)
    else:
        parser.print_help()

if __name__ == "__main__":
    main() 