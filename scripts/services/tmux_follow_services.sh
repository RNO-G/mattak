#!/bin/bash
#
# tmux_follow_services.sh
#
# Follow the journalctl logs of the rno-g-autoconverter systemd services
# for all stations in a tmux session named "autoconverter".
#
# One pane is created per station, each running
#   journalctl -u rno-g-autoconverter@<station> -f
#
# By default all panes are tiled in a single window. With -p N the panes
# are distributed over multiple windows ("pages") with at most N panes
# each; new windows are created as needed and named after the station
# ids they contain.
#
# If the session already exists, the script simply attaches to it.
# When any pane dies, the whole session is killed.
#
# Usage: tmux_follow_services.sh [-p PANES_PER_WINDOW] [-h]

STATIONS="11 12 13 14 21 22 23 24 15 25 34 35"
SESSION="autoconverter"

# usage: Print help text describing the script and its options.
usage() {
  cat <<EOF
Usage: $(basename "$0") [-p PANES_PER_WINDOW] [-h]

Follow the journalctl logs of the rno-g-autoconverter services for all
stations ($STATIONS) in a tmux session named "$SESSION".

Options:
  -p PANES_PER_WINDOW  Maximum number of panes per tmux window ("page").
                       New windows are created as needed and named after
                       the station ids they contain.
                       Default: 0 (all panes in a single window).
  -h                   Show this help and exit.
EOF
}

# setup_window WIN: Configure window WIN to show pane titles in the
# pane borders.
setup_window() {
  tmux set-option -w -t "$1" pane-border-status top
  tmux set-option -w -t "$1" pane-border-format "#{pane_title} "
}

# start_follow STATION WIN: Title the active pane of window WIN and start
# following the service log of STATION in it.
start_follow() {
  tmux select-pane -t "$2" -T "Station $1"
  tmux send-keys -t "$2" "journalctl -u rno-g-autoconverter@$1 -f" Enter
}

# --- Option parsing ---------------------------------------------------------

panes_per_window=0  # 0 means: no limit, all panes in a single window

while getopts ":p:h" opt; do
  case $opt in
    p) panes_per_window=$OPTARG ;;
    h) usage; exit 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument" >&2; usage >&2; exit 1 ;;
  esac
done

if ! [[ "$panes_per_window" =~ ^[0-9]+$ ]]; then
  echo "-p expects a non-negative integer, got '$panes_per_window'" >&2
  exit 1
fi

# --- Sanity checks ----------------------------------------------------------

# Don't nest tmux sessions
if [ -n "$TMUX" ]; then
  echo "Already inside a tmux session, refusing to nest"
  exit 1
fi

# Don't recreate session if it already exists
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "Session $SESSION already exists, attaching..."
  tmux attach -t "$SESSION"
  exit 0
fi

# --- Session setup ----------------------------------------------------------

# Create the session and remember the id of its initial window
win=$(tmux new-session -d -s "$SESSION" -P -F "#{window_id}")
first_win=$win
setup_window "$win"

panes_in_window=0  # panes created so far in the current window
win_stations=""    # station ids shown in the current window (for its name)

for st in $STATIONS; do
  if [ "$panes_in_window" -eq 0 ]; then
    :  # first pane of the first window already exists
  elif [ "$panes_per_window" -gt 0 ] && [ "$panes_in_window" -ge "$panes_per_window" ]; then
    # Current window is full, open a new one
    win=$(tmux new-window -t "$SESSION" -P -F "#{window_id}")
    setup_window "$win"
    panes_in_window=0
    win_stations=""
  else
    tmux split-window -t "$win"
  fi

  start_follow "$st" "$win"
  tmux select-layout -t "$win" tiled

  panes_in_window=$((panes_in_window + 1))
  win_stations="${win_stations:+$win_stations,}$st"

  # With paging enabled, name each window after the stations it contains
  if [ "$panes_per_window" -gt 0 ]; then
    tmux rename-window -t "$win" "St $win_stations"
  fi
done

# Kill the whole session as soon as any pane dies
tmux set-hook -t "$SESSION" pane-died "kill-session -t $SESSION"

tmux select-window -t "$first_win"
tmux attach -t "$SESSION"
