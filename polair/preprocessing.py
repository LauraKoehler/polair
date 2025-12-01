#import _helpers as h
import xarray as xr
import pandas as pd

def configure_preprocessing_parser(parser):
    parser.add_argument(
        "-f",
        "--flight",
        help="Research flight number (integer)",
        default=None,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--config",
        metavar="CONFIG_FILE",
        help="configuration file (yaml)",
        default=None,
        required=True,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        metavar="DEBUG",
        help='Set the level of verbosity [DEBUG, INFO," " WARNING, ERROR]',
        required=False,
        default="INFO",
    )

    parser.set_defaults(func=run)
    
def run(args):
    flight = args.flight
    config_file = args.config
    print(flight)