from typing import Iterator, List

class TabixError(Exception):
    ...

class TabixFile:
    def query(self, chrom: str, start: int, end: int) -> Iterator[List[str]]:
        ...

def open(filename: str) -> TabixFile:
    ...
