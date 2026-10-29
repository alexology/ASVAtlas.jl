#!/bin/sh
# Writes the usage reports served (with a password) at https://<domain>/stats/:
#   stats/index.html  visits, from the Caddy access log (GoAccess); the history
#                     is kept in stats/db, so rotated logs are not needed
#   stats/usage.txt   searches and networks (ASVAtlas usage log, counts only)
# Run it in /srv/asvatlas, e.g. every night with cron (crontab -e):
#   15 3 * * * cd /srv/asvatlas && sh update-stats.sh > /dev/null 2>&1
set -e
cd "$(dirname "$0")"
mkdir -p stats/db

if [ -s logs/access.log ]; then
    docker run --rm -v "$PWD/logs:/logs:ro" -v "$PWD/stats:/stats" allinurl/goaccess \
        /logs/access.log --log-format=CADDY --anonymize-ip --no-query-string \
        --persist --restore --db-path=/stats/db -o /stats/index.html
fi

docker compose exec -T asvatlas julia /opt/ASVAtlas.jl/bin/asvatlas-usage.jl /opt/stats/usage.tsv \
    > stats/usage.txt 2> /dev/null || echo "No search yet." > stats/usage.txt
