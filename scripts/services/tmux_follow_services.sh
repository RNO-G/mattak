#!/bin/bash
STATIONS="23 11 22 12 21 13 24 14"
SESSION="autoconverter"

# Don't nest tmux sessions
if [ -n "$TMUX" ]; then
  echo "Already inside a tmux session, refusing to nest"
  exit 1
fi

# Don't recreate session if it already exists
if tmux has-session -t $SESSION 2>/dev/null; then
  echo "Session $SESSION already exists, attaching..."
  tmux attach -t $SESSION
  exit 0
fi


tmux new-session -d -s $SESSION

# Show pane titles as borders
tmux set-option -t $SESSION pane-border-status top
tmux set-option -t $SESSION pane-border-format " Station #{pane_title} "

i=0
for st in $STATIONS; do
  if [ $i -eq 0 ]; then
    tmux select-pane -t $SESSION -T "Station ${st}"
    tmux send-keys -t $SESSION "journalctl -u rno-g-autoconverter@${st} -f" Enter
  else
    tmux split-window -t $SESSION
    tmux select-pane -t $SESSION -T "Station ${st}"
    tmux send-keys -t $SESSION "journalctl -u rno-g-autoconverter@${st} -f" Enter
  fi
  tmux select-layout -t $SESSION tiled
  i=$((i+1))
done


tmux set-hook -t $SESSION pane-died "kill-session -t $SESSION"

tmux attach -t $SESSION
