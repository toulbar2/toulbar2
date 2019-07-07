#!/bin/bash
git pull origin kad 1>/dev/null 2>&1
git add . 1>/dev/null 2>&1
git commit -m "master-worker" 1>/dev/null 2>&1
git push origin kad 1>/dev/null 2>&1
git push origin kad 1>/dev/null 2>&1
git push origin kad | grep "Everything up-to-date"
sleep 1s
