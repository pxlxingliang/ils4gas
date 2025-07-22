import argparse
import os

from ils4gas.env import set_envs, create_workpath
from ils4gas.google_adk_template import gen_google_adk_template

from importlib.metadata import version
__version__ = version("ils4gas")

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="AbacusAgent Command Line Interface")
    
    parser.add_argument(
        "--transport",
        type=str,
        default=None,
        choices=["sse", "streamable-http"],
        help="Transport protocol to use (default: sse), choices: sse, streamable-http"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=None,
        help="Port to run the MCP server on (default: 50001)"
    )
    parser.add_argument(
        "--host",
        type=str,
        default=None,
        help="Host to run the MCP server on (default: localhost)"
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the version of Ils4gasAgentTools"
    )
    parser.add_argument(
        "--template",
        action="store_true",
        help="Generate a Google ADK template for ils4gas"
    )
    
    args = parser.parse_args()
    
    return args

def print_version():
    """
    Print the version.
    """
    print(f"\nIls4gasAgentTools Version: {__version__}")
    print("For more information, visit: https://gitee.com/pxlxingliang/ils4gas\n")
    
def print_address():
    """
    Print the address of the MCP server based on environment variables.
    """
    address = f"{os.environ['ILS4GAS_HOST']}:{os.environ['ILS4GAS_PORT']}"
    if os.environ["ILS4GAS_TRANSPORT"] == "sse":
        print("Address:", address + "/sse")
    elif os.environ["ILS4GAS_TRANSPORT"] == "streamable-http":
        print("Address:", address + "/mcp")
    else:
        raise ValueError("Invalid transport protocol specified. Use 'sse' or 'streamable-http'.")

def main():
    """
    Main function to run the MCP tool.
    """
    print_version()
    args = parse_args()  

    set_envs(
        transport_input=args.transport,
        port_input=args.port, 
        host_input=args.host)
    create_workpath()
    
    if args.version:
        print_version()
        return
    
    if args.template:
        gen_google_adk_template()
        return
    
    from ils4gas.init_mcp import mcp
    from ils4gas.env import load_tools
    load_tools()  

    print_address()
    mcp.run(transport=os.environ["ILS4GAS_TRANSPORT"])

if __name__ == "__main__":
    main()
