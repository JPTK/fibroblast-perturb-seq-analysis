from pathlib import Path

DIR_DATA_RNA = Path("/data/torsten/GSE261783_RAW/")
DIR_DATA_GENERATED = Path("/data/torsten/data_generated/")
DIR_PLOTS = Path("/data/torsten/plots/")

for dir in [DIR_DATA_RNA, DIR_DATA_GENERATED, DIR_PLOTS]:
    dir.mkdir(exist_ok=True, parents=True)