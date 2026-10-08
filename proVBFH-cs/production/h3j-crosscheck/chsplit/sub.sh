#!/bin/bash
# sbatch with the chsplit defaults: nice 50, bad nodes and et07-40 excluded
EXCL="$(paste -sd, /ptmp/mpp/akarlber/cs-production/bad_nodes),et[07-40]"
exec sbatch --nice=50 --exclude="$EXCL" "$@"
