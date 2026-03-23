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
