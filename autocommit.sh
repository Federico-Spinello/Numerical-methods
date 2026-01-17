echo "CRON RUN $(date)" >> /home/federico/cron-test.log
#!/bin/bash

cd "$(dirname "$0")" || exit 1

git add -A

if git diff --cached --quiet; then
    exit 0
fi

git commit -m "auto: $(date '+%Y-%m-%d %H:%M:%S')"
git push

