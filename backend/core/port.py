import socket


def find_free_port(port: int, host: str = "0.0.0.0", max_attempts: int = 100) -> int:
    for offset in range(max_attempts):
        p = port + offset
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((host, p))
                return p
            except OSError:
                continue
    raise RuntimeError(f"No free port found starting from {port}")
