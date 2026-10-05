#!/bin/bash
# Hourly: status, and feed the remaining hxswg136 NLO exclusive jobs.
P=/ptmp/mpp/akarlber/cs-production; Q=$P/prod-lonlo; S=$(dirname "$0")
"$S/status.sh"
# NNLO exclusive: keep at most ~2,500 of my jobs in the queue (AK, 5 Oct:
# 40k jobs queueing is not sustainable); never-started lines are fed in
for s in hxswg136 p1506; do
    "$S/feed_missing.sh" $s-excl-nnlo $P/bin/proVBFH-cs-$s $P/prod-nnlo/$s/excl/jobs.list 12:00:00 2500 500
done
"$S/feed.sh" hxswg136-nlo-excl $P/bin/proVBFH-cs-hxswg136 $Q/hxswg136/nlo-excl/jobs.list 6:00:00
# resubmit unfinished tasks of lists whose arrays have ended
for s in p1506 hxswg136; do
    "$S/resubmit.sh" $s-incl-nnlo $P/bin/proVBFH-cs-$s $P/prod-nnlo/$s/incl/jobs.list 2:00:00
    "$S/resubmit.sh" $s-nlo-incl $P/bin/proVBFH-cs-$s $Q/$s/nlo-incl/jobs.list 6:00:00
    "$S/resubmit.sh" $s-lo-incl $P/bin/proVBFH-cs-$s $Q/$s/lo-incl/jobs.list 6:00:00
done
"$S/resubmit.sh" p1506-nlo-excl $P/bin/proVBFH-cs-p1506 $Q/p1506/nlo-excl/jobs.list 6:00:00
[ "$(cat $Q/hxswg136/nlo-excl/jobs.list.next)" -gt 6600 ] && \
    "$S/resubmit.sh" hxswg136-nlo-excl $P/bin/proVBFH-cs-hxswg136 $Q/hxswg136/nlo-excl/jobs.list 6:00:00
true
