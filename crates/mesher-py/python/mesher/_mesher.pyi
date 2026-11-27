"""Type stubs for the extension module (``src/python.rs``)."""

from os import PathLike
from typing import Callable, Dict, List, Optional, Tuple, Union

__version__: str

def _run(
    job_json: str,
    worker: Union[str, PathLike[str]],
    timeout: float,
    memory: Optional[int] = None,
    progress: Optional[Callable[[str, float], None]] = None,
) -> Tuple[str, Optional[str]]: ...
def _run_request(
    request_json: str,
    worker: Union[str, PathLike[str]],
    timeout: float,
    memory: Optional[int] = None,
    progress: Optional[Callable[[str, float], None]] = None,
) -> Tuple[str, Optional[str]]: ...
def _build_job(request_json: str) -> Optional[str]: ...
def _next_job(report_json: str, output_dir: str, min: str) -> Optional[Tuple[str, int, int]]: ...
def _read_mesh(path: Union[str, PathLike[str]]) -> Dict[str, Union[bytes, List[str]]]: ...
