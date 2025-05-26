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
import time # For timing
from vcf_agent import metrics # Import metrics module

# OpenTelemetry Imports
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider as SdkTracerProvider
from .tracing import init_tracer, setup_auto_instrumentation

# Initialize OpenTelemetry Tracer for the CLI
# This should be done once, as early as possible.
# The service name helps distinguish traces from different components.
cli_tracer = init_tracer(service_name="vcf-agent-cli")
setup_auto_instrumentation() # Placeholder, will configure auto-instrumentation

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
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")
    subparsers.required = False # Allow no subcommand to print help

    # Common arguments for agent interaction (can be added to main or specific subcommands)
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
    parser.add_argument(
        "--ollama-model", type=str, default=None, # Default is None, SessionConfig has its own default
        help="Specify the Ollama model name (e.g., qwen:4b, llama2). Overrides SessionConfig default."
    )

    # ADD 'ask' SUBCOMMAND
    ask_parser = subparsers.add_parser("ask", help="Ask the VCF agent a question or give it a command. This is the default if no other subcommand is given and text follows.")
    ask_parser.add_argument("prompt_text", type=str, help="The prompt text for the agent.")

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

    # Kuzu: populate from VCF
    parser_populate_kuzu = subparsers.add_parser("populate-kuzu-from-vcf", help="Populate Kuzu graph database from a VCF file.")
    parser_populate_kuzu.add_argument("--vcf_file_path", type=str, required=True, help="Path to the VCF file.")
    parser_populate_kuzu.add_argument("--kuzu_db_path", type=str, default="./kuzu_db", help="Path to Kuzu database directory.")
    # Optionally, add --sample_name_override if needed for the CLI too

    args = parser.parse_args()

    # If no command is given, but there are unparsed args, assume it's an 'ask' command
    # This makes 'ask' the default if non-option arguments are present.
    # We need to re-parse if this happens.
    if args.command is None and any(arg for arg in sys.argv[1:] if not arg.startswith('-')):
        # Check if the first non-option argument is a known command, if so, user forgot to specify it.
        # This is a basic check. More robust would be to inspect sys.argv directly.
        potential_command = next((arg for arg in sys.argv[1:] if not arg.startswith('-')), None)
        if potential_command not in subparsers.choices:
            # Re-parse with 'ask' prepended if it's not another known command
            # This assumes all non-flag args constitute the prompt
            prompt_args = [arg for arg in sys.argv[1:] if not arg.startswith('--')]
            other_args = [arg for arg in sys.argv[1:] if arg.startswith('--')]
            
            # Filter out any already parsed known global options from prompt_args
            # This is tricky because argparse has already consumed them.
            # A simpler approach for now: if no command, and sys.argv[1] isn't an option, it's a prompt.
            
            # Let's try a simpler default handling for now without re-parsing.
            # If args.command is None, and there was something that could be a prompt,
            # we'll manually set args.command to 'ask' and construct prompt_text.
            # This requires sys.argv inspection before parse_args, or a try-catch parse.

            # The following logic is based on the idea that if 'command' is None
            # AND there were positional arguments supplied that weren't consumed by options,
            # these positional arguments are intended for the 'ask' command.
            # This is hard to do robustly AFTER parse_args().
            # The `subparsers.required = False` and a check later is better.
            pass # Will handle default action later

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
    elif args.command == "populate-kuzu-from-vcf":
        with cli_tracer.start_as_current_span("cli.populate_kuzu_from_vcf") as populate_span:
            populate_span.set_attribute("vcf.file_path", args.vcf_file_path)
            populate_span.set_attribute("kuzu.db_path", args.kuzu_db_path)

            from vcf_agent.graph_integration import get_kuzu_db_connection, create_schema
            from vcf_agent.vcf_utils import populate_kuzu_from_vcf
            import gc

            command_name = "populate-kuzu-from-vcf"
            status = "success"
            start_time = time.time()
            error_observed = None
            # Ensure CLI metrics are defined; they are not auto-registered to http_registry in metrics.py
            # For push, we collect their current state.
            cli_duration = metrics.CLI_COMMAND_DURATION_SECONDS
            cli_requests = metrics.CLI_COMMAND_REQUESTS_TOTAL
            cli_errors = metrics.CLI_COMMAND_ERRORS_TOTAL

            try:
                kuzu_conn = None # Define kuzu_conn here to ensure it's in scope for finally
                try:
                    kuzu_conn = get_kuzu_db_connection(db_path=args.kuzu_db_path)
                    create_schema(kuzu_conn) # Ensure schema exists
                    
                    metrics.log.info(f"Populating Kuzu graph at '{args.kuzu_db_path}' from VCF file '{args.vcf_file_path}'...")
                    counts = populate_kuzu_from_vcf(kuzu_conn=kuzu_conn, vcf_path=args.vcf_file_path)
                    metrics.log.info(f"Successfully populated Kuzu graph from VCF. Variants added: {counts.get('variants', 0)}, Samples added/found: {counts.get('samples', 0)}, Links created: {counts.get('links', 0)}. ")
                    populate_span.set_attribute("kuzu.variants_added", counts.get('variants', 0))
                    populate_span.set_attribute("kuzu.samples_added_found", counts.get('samples', 0))
                    populate_span.set_attribute("kuzu.links_created", counts.get('links', 0))

                except FileNotFoundError as e:
                    status = "error"
                    error_observed = e
                    metrics.log.error(f"Error: VCF file not found at '{args.vcf_file_path}'. Details: {e}", file=sys.stderr)
                    populate_span.record_exception(e)
                    populate_span.set_status(trace.StatusCode.ERROR, f"VCF file not found: {args.vcf_file_path}")
                except Exception as e:
                    status = "error"
                    error_observed = e
                    metrics.log.error(f"Error during Kuzu population: {e}", file=sys.stderr)
                    populate_span.record_exception(e)
                    populate_span.set_status(trace.StatusCode.ERROR, "Error during Kuzu population")
                finally:
                    if kuzu_conn:
                        del kuzu_conn 
                    gc.collect() 
            finally:
                duration = time.time() - start_time
                cli_duration.labels(command=command_name, status=status).observe(duration)
                cli_requests.labels(command=command_name, status=status).inc()
                
                metrics_to_push_list = [cli_duration, cli_requests]

                if status == "error" and error_observed:
                    error_type = type(error_observed).__name__
                    cli_errors.labels(command=command_name, error_type=error_type).inc()
                    metrics_to_push_list.append(cli_errors)
                
                metrics.push_job_metrics(job_name=command_name, metrics_to_push=metrics_to_push_list)

                if error_observed: # Re-raise or exit after pushing metrics
                    sys.exit(1)
        return

    # Agent interaction logic (now specifically for 'ask' or default)
    # This block should be executed if args.command == "ask"
    # OR if args.command is None (meaning no other subcommand was explicitly called)
    
    prompt_to_run = None
    if args.command == "ask":
        prompt_to_run = args.prompt_text
    elif args.command is None:
        # This is the tricky part: how to get the prompt if no subcommand was specified?
        # parse_args() with nargs='?' for a positional argument on main parser
        # normally handles this.
        # For now, if command is None, we assume the user wanted help or made an error.
        # The 'ask' command must be explicit.
        parser.print_help()
        sys.exit(0) # Exit cleanly after help

    # This means args.command MUST be "ask" to reach here
    # (or other future agent-related subcommands)
    
    if prompt_to_run: # This will be true if command is 'ask'
        # --- MINIMAL OTEL TEST WAS HERE ---
        # print("[CLI] Starting MINIMAL OTEL TEST.")
        # with cli_tracer.start_as_current_span("cli.minimal_otel_test_span") as test_span:
        #     test_span.set_attribute("test.attribute", "hello_jaeger")
        #     print(f"[CLI] Minimal test span created: trace_id={test_span.get_span_context().trace_id:x}, span_id={test_span.get_span_context().span_id:x}")
        #     time.sleep(1) # Give a moment for things to happen
        # print("[CLI] Minimal test span ended.")
        # --- END MINIMAL OTEL TEST ---

        # UNCOMMENT THE REST OF THE 'ask' LOGIC FOR THIS TEST
        # """ # Remove this line to uncomment
        if args.raw:
            os.environ["VCF_AGENT_RAW_MODE"] = "1"

        from vcf_agent.agent import get_agent_with_session
        from vcf_agent.config import SessionConfig
        
        model_provider = cast(Literal["ollama", "openai", "cerebras"], args.model)
        session_config = SessionConfig(
            raw_mode=args.raw if args.raw else None,
            model_provider=model_provider,
            credentials_file=args.credentials,
            reference_fasta=args.reference,
            ollama_model_name=args.ollama_model if args.ollama_model else None # Pass it here
        )
        
        agent_instance = get_agent_with_session( # Renamed to avoid conflict
            session_config=session_config,
            model_provider=model_provider
        )
        
        # The actual agent call within a span
        with cli_tracer.start_as_current_span("cli.agent_ask_interaction") as prompt_span: # Renamed span
            prompt_span.set_attribute("agent.prompt", prompt_to_run)
            prompt_span.set_attribute("agent.model_provider", model_provider)
            prompt_span.set_attribute("agent.raw_mode", args.raw if args.raw else False)
            
            current_span_obj_before_agent = trace.get_current_span()
            current_ctx_before_agent = current_span_obj_before_agent.get_span_context()
            if current_ctx_before_agent and current_ctx_before_agent.is_valid:
                print(f"[CLI] Before agent call. Current span: trace_id={hex(current_ctx_before_agent.trace_id)}, span_id={hex(current_ctx_before_agent.span_id)}")
            else:
                print("[CLI] Before agent call. No valid current span or context.")

            try:
                raw_agent_response = agent_instance(prompt_to_run)
                
                output_text = ""
                if isinstance(raw_agent_response, str):
                    output_text = raw_agent_response
                elif hasattr(raw_agent_response, 'response') and isinstance(getattr(raw_agent_response, 'response'), str):
                    output_text = getattr(raw_agent_response, 'response')
                elif hasattr(raw_agent_response, 'text') and isinstance(getattr(raw_agent_response, 'text'), str):
                    # Common alternative attribute name for response text
                    output_text = getattr(raw_agent_response, 'text')
                else:
                    output_text = str(raw_agent_response) 
                
                prompt_span.set_attribute("agent.response_length", len(output_text) if output_text is not None else 0)
                prompt_span.set_status(trace.StatusCode.OK)
                print(output_text)
            except Exception as e:
                prompt_span.record_exception(e)
                prompt_span.set_status(trace.StatusCode.ERROR, str(e))
                print(f"Error during agent execution: {e}", file=sys.stderr)
                sys.exit(1) # Exit with error after logging
        
        # Metrics server start - consider if this should only be for 'ask'
        metrics.start_metrics_http_server() 
        # """ # Remove this line to uncomment
        # After the command has executed, force flush any pending spans
        current_provider = trace.get_tracer_provider()
        if isinstance(current_provider, SdkTracerProvider):
            print("[CLI] Forcing flush of OTel spans...")
            try:
                current_provider.force_flush(timeout_millis=5000) # 5 second timeout
                print("[CLI] OTel spans flushed.")
            except Exception as e:
                print(f"[CLI] Error during OTel span flush: {e}")
        else:
            print("[CLI] No SDK TracerProvider found, skipping flush.")
            
    else:
        parser.print_help()

if __name__ == "__main__":
    main() 