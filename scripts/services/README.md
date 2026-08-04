These file must be installed, everything is set up in the Makefile

`sudo make install-systemd`

then do
`sudo systemctl restart rno-g-autoconverter.target`

use `tmux_follow_services.sh`to follow the logs of all stations together in a tmux session.

Usefull commands:

Install 
`sudo systemctl daemon-reload`

Start everything
`sudo systemctl start rno-g-autoconverter.target`

Watch a specific station (e.g 23)
`journalctl -u rno-g-autoconverter@23 -f`

Check status of all stations 
`systemctl status 'rno-g-autoconverter@*'`

Stop everything  
`systemctl stop rno-g-autoconverter.target`

Restart a single station 
`sudo systemctl restart rno-g-autoconverter@23` 

See logs for all stations together 
`journalctl -u 'rno-g-autoconverter@*' -f`

## Rootified cleanup

Purges the bulk `*.root` files, `cfg/`, and logs of rootified runs
(`/data/rootified/station*/run*`) older than 30 days (by last modification),
once a day. `rootified-list.txt` and `aux/acq-file-list.txt` are
deliberately kept: they are exactly what `rno-g-convert-run` checks to
decide a run is already up to date, so keeping them stops the autoconverter
from reconverting a purged run from scratch. `aux/comment.txt` and
`aux/runinfo.txt` are kept too (tiny, useful metadata); everything else in
`aux/` (e.g. `flower_end.json.gz`, which can be large) is purged along with
the rest. To protect a run from being purged at all, create a `.keep` file
inside its run directory, e.g.:
`touch /data/rootified/station23/run1144/.keep`

Install
`sudo make install-systemd-cleanup`

Enable and start the daily timer
`sudo systemctl enable --now rno-g-cleanup-rootified.timer`

Check when it last/next runs
`systemctl status rno-g-cleanup-rootified.timer`

See logs
`journalctl -u rno-g-cleanup-rootified -f`

Run it manually (dry-run, doesn't delete anything)
`./scripts/daq/rno-g-cleanup-rootified 30 /data/rootified 1`
