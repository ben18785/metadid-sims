# Metapsy datasets

## Source

Metapsy project's public repositories: `github.com/metapsy-project/data-<slug>`.

## What is in there

2,699 rows, metadata columns first, then data.

| design | rows |
| --- | --- |
| `did` — pre and post, both arms | 2,075 |
| `rct` — post only, both arms | 609 |
| `did_change` — change scores only | 15 |

`grief`, `panic` and `psychosis` are entirely post-only.

## Limits — read before fitting

**Some rows are not independent.** 20% share a control arm with another row in the
same study: a three-arm trial contributes each treatment against the same
control, with the control numbers repeated verbatim. Pass `study_id` as the
study, not `Intervention_ID`, or the shared controls are double-counted and
uncertainty is understated.

**Timepoint is not carried.** Metapsy distinguishes `post` from `FU1`/`FU2`;
so a study measured at both appears twice and a `did` row may be either.
Where post-test cells are empty the two rows become byte-identical; those are
among the dropped rows, but the ambiguity remains for the rows that survive.

**Some studies have baselines and no outcome.** Metapsy stored a Hedges' *g*
derived from some other statistic (`outcome_type = "other stat"`) and never
held raw post-test means. Such rows are dropped.

**661 randomisation values are blank.** Noting that `metadid` will in these
cases infer randomisation as `none` or not randomised with the default setting.
Future work could go back to the sources and add randomisation data which was
not available in Metapsy.