import socket

from peroxide_engine import cli


def test_free_port_steps_past_a_port_in_use():
    with socket.socket() as held:
        held.bind(("", 0))
        held.listen()
        taken = held.getsockname()[1]
        assert cli.free_port(taken, tries=20) > taken
