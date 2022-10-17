#!/usr/bin/env python3
from pathlib import Path
import argparse

import acts
import acts.examples


u = acts.UnitConstants


<<<<<<< Updated upstream
def RootBFieldWrite(bField, fileName, treeName="solenoid", level=acts.logging.VERBOSE):
    cfg = acts.examples.RootBFieldWriter.Config()
    cfg.bField = bField
    cfg.gridType = acts.examples.RootBFieldWriter.GridType.rz
    cfg.fileName = str(fileName)
    cfg.treeName = treeName
    acts.examples.RootBFieldWriter.run(cfg, level)
    return cfg


def CsvBFieldWrite(bField, fileName, level=acts.logging.VERBOSE):
    cfg = acts.examples.CsvBFieldWriter.ConfigXyzGridless()
    cfg.bField = bField
    cfg.range = ((-10000., 10000.), (-10000., 10000.), (-15000., 15000.))
    cfg.bins = (201, 201, 301)
    cfg.fileName = str(fileName)
    acts.examples.CsvBFieldWriter.runXyzGridless(cfg, level)
    return cfg


def runBFieldWriting(outputDir: Path, rewrites: int = 0):
=======
def runRootBFieldWriting(args: argparse.Namespace):
    if args.useXyz:
        raise ValueError("ROOT B-field writer only supports cylindrical coordinates.")

>>>>>>> Stashed changes
    solenoid = acts.SolenoidBField(
        radius=1200 * u.mm, length=6000 * u.mm, bMagCenter=2 * u.T, nCoils=1194
    )

    if args.bins is None:
        bins = (10, 10)
    else:
        bins = (args.bins[0] or 10, args.bins[1] or 10)

    rlim = (
        ((args.min[0] or 0) if args.min else 0) * u.mm,
        ((args.max[0] or 1200.0) if args.max else 1200) * u.mm,
    )

    zlim = (
        ((args.min[1] or -5000.0) if args.min else -5000.0) * u.mm,
        ((args.max[1] or 5000.0) if args.max else 5000.0) * u.mm,
    )

    field = acts.solenoidFieldMap(
        rlim=rlim,
        zlim=zlim,
        nbins=bins,
        field=solenoid,
    )

    print("Solenoid ready")

<<<<<<< Updated upstream
    cfg = RootBFieldWrite(field, outputDir / "solenoid.root")
    CsvBFieldWrite(field, outputDir / "solenoid.csv")
=======
    cfg = acts.examples.RootBFieldWriter.Config()
    cfg.bField = field
    cfg.gridType = acts.examples.RootBFieldWriter.GridType.rz
    cfg.fileName = str(args.output / "solenoid.root")
    cfg.treeName = "solenoid"

    acts.examples.RootBFieldWriter.run(cfg, acts.logging.VERBOSE)
>>>>>>> Stashed changes

    for i in range(args.rewrites):
        print(f"Now read back {cfg.fileName}")

        field2 = acts.examples.MagneticFieldMapRz(cfg.fileName, tree="solenoid")
<<<<<<< Updated upstream
        cfg2 = RootBFieldWrite(field2, outputDir / f"solenoid{i+2}.root")
        CsvBFieldWrite(field2, outputDir / f"solenoid{i+2}.csv")
=======

        cfg2 = acts.examples.RootBFieldWriter.Config()
        cfg2.bField = field2
        cfg2.gridType = acts.examples.RootBFieldWriter.GridType.rz
        cfg2.fileName = str(args.output / f"solenoid{i+2}.root")
        cfg2.treeName = "solenoid"

        acts.examples.RootBFieldWriter.run(cfg2, acts.logging.VERBOSE)

>>>>>>> Stashed changes
        cfg = cfg2

    print("Done")


def runCsvBFieldWriting(args: argparse.Namespace):
    if args.bins is None or any(x is None for x in args.bins):
        raise ValueError("All bins must be specified for CSV writing.")

    solenoid = acts.SolenoidBField(
        radius=1200 * u.mm, length=6000 * u.mm, bMagCenter=2 * u.T, nCoils=1194
    )

    if args.useXyz:
        cfg = acts.examples.CsvBFieldWriter.ConfigXyzGridless()
        cfg.bField = solenoid
        cfg.fileName = str(args.output)
        cfg.bins = args.bins
        cfg.range = ((-10000, 10000), (-10000, 10000), (-15000, 15000))

        acts.examples.CsvBFieldWriter.runXyzGridless(cfg, acts.logging.VERBOSE)
    else:
        cfg = acts.examples.CsvBFieldWriter.ConfigRzGridless()
        cfg.bField = solenoid
        cfg.fileName = str(args.output)
        cfg.bins = args.bins
        cfg.range = ((0, 10000), (-15000, 15000))

        acts.examples.CsvBFieldWriter.runRzGridless(cfg, acts.logging.VERBOSE)


class RangeAction(argparse.Action):
    def __init__(self, **kwargs):
        range_type = kwargs.pop("range_type", int)

        super(RangeAction, self).__init__(**kwargs)

        self._range_type = range_type

    def __call__(self, parser, namespace, values, option_string=None):
        if not hasattr(namespace, "useXyz"):
            raise ValueError("Coordinate system must be set before range.")

        num = 3 if namespace.useXyz else 2

        vals = values.split(":")

        if len(vals) != num:
            raise ValueError("Range must have exactly " + str(num) + " values.")

        setattr(
            namespace,
            self.dest,
            [self._range_type(v) if v != "" else None for v in vals],
        )


if "__main__" == __name__:
    parser = argparse.ArgumentParser()

    coordinate_group = parser.add_mutually_exclusive_group(required=True)
    coordinate_group.add_argument(
        "--xyz",
        dest="useXyz",
        action="store_true",
        help="use Cartesian coordinate system",
    )
    coordinate_group.add_argument(
        "--rz",
        dest="useXyz",
        action="store_false",
        help="use cylindrical coordinate system",
    )

    parser.add_argument(
        "-b",
        "--bins",
        action=RangeAction,
        range_type=int,
        help="number of bins in result",
        metavar="[a1]:...:[an]",
    )
    parser.add_argument(
        "--min",
        action=RangeAction,
        range_type=float,
        help="minimum coordinate",
        metavar="[c1]:...:[cn]",
    )
    parser.add_argument(
        "--max",
        action=RangeAction,
        range_type=float,
        help="maximum coordinate",
        metavar="[c1]:...:[cn]",
    )

    subparsers = parser.add_subparsers(
        help="format to use for output",
        required=True,
    )

    parser_root = subparsers.add_parser("root", help="lol")
    parser_root.add_argument(
        "-r",
        "--rewrites",
        type=int,
        default=0,
        metavar="N",
        help="number of times to rewrite file",
    )
    parser_root.add_argument(
        "--output",
        type=Path,
        default=Path.cwd(),
        help="output directory to write to",
    )
    parser_root.set_defaults(func=runRootBFieldWriting)

    parser_csv = subparsers.add_parser("csv")
    parser_csv.add_argument(
        "--output",
        type=Path,
        default="bfield.csv",
        help="output file to write to",
    )
    parser_csv.set_defaults(func=runCsvBFieldWriting)

    args = parser.parse_args()
    print(args)
    args.func(args)
