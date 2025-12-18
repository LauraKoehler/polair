import logging

from .preprocessing import configure_preprocessing_parser
from .noseboom import configure_noseboom_parser
from ._version import __version__


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-v",
        "--verbose",
        metavar="DEBUG",
        help="Set the level of verbosity [DEBUG, INFO, WARNING, ERROR]",
        required=False,
        default="INFO",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=__version__,
    )

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers(dest="subcommand")
    configure_preprocessing_parser(subparsers.add_parser("preprocessing"))
    configure_noseboom_parser(subparsers.add_parser("noseboom"))

    return parser


def main():
    args = get_parser().parse_args()

    logging.basicConfig(level=logging.getLevelName(args.verbose))

    return args.func(args)


if __name__ == "__main__":
    exit(main())
