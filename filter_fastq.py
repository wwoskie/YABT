import argparse

from yet_another_bioinformatic_tool import FastQFiltrator


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter fastq files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_path",
        "-i",
        help="Path to input",
        type=str,
    )

    parser.add_argument(
        "--output_path",
        "-o",
        help="Path to output",
        default=None,
        type=str,
    )

    parser.add_argument(
        "--gc_lower_bound",
        "-gcl",
        help="GC content lower bound",
        default=0,
        type=float,
    )

    parser.add_argument(
        "--gc_upper_bound",
        "-gcu",
        help="GC content upper bound",
        default=1,
        type=float,
    )

    parser.add_argument(
        "--length_lower_bound",
        "-ll",
        help="Length lower bound",
        default=0,
        type=int,
    )

    parser.add_argument(
        "--length_upper_bound",
        "-lu",
        help="Length upper bound",
        default=2**32,
        type=int,
    )

    parser.add_argument(
        "--quality_threshold",
        "-qt",
        help="Length upper bound",
        default=0,
        type=float,
    )

    parser.add_argument(
        "--logs_dir",
        "-ld",
        help="Path to directory for logs",
        default="",
        type=str,
    )

    parser.add_argument(
        "--rewrite",
        "-r",
        help="If should rewrite output",
        action="store_true",
    )
    parser.set_defaults(rewrite=False)

    args = parser.parse_args()

    filtrator = FastQFiltrator(
        path_to_input=args.input_path,
        gc_bounds=(args.gc_lower_bound, args.gc_upper_bound),
        length_bounds=(args.length_lower_bound, args.length_upper_bound),
        quality_threshold=args.quality_threshold,
        logs_dir=args.logs_dir,
    )

    filtrator.filter_fastq()
    filtrator.write_to_file(args.output_path, rewrite=args.rewrite)


if __name__ == "__main__":
    main()
