# Understanding your results

### The metaCpG score

The **metaCpG score** is a single number summarising the methylation pattern of
your sample across a small panel of prognostic CpG sites. It is calculated as a
weighted average of the methylation signal at those sites and then standardised
against our discovery cohort, so that a score of **0** is the mean, **positive**
scores indicate a more adverse methylation profile, and **negative** scores a
more favourable one. The distribution plot shows where your sample falls relative
to all previously analysed cases.

### Predicted overall survival

The risk score is combined with the relevant clinical features and entered into a
survival model to produce a **predicted overall survival (OS)** — the estimated
probability that a patient with this molecular and clinical profile would be
alive at a given time point.

- **SHH tumours** — a predicted **5-year OS**, adjusted for *MYCN* amplification
  status.
- **Group 3/4 tumours** — a predicted **5-year OS**, adjusted for metastatic
  (M) stage and stratified by *MYC* amplification status.

### The late metaCpG (Group 3/4)

For Group 3/4 tumours, we additionally report a **late metaCpG** estimate: a
predicted probability of a **late event (death by 10 years)**. Some Group 3/4
tumours that do well in the first few years remain at risk of later relapse, and
this second signature is designed to flag that longer-term risk that the early
5-year estimate may not capture. The two readouts are complementary — read them
together rather than in isolation.

### How to read the numbers

A predicted 5-year OS of, say, **80%** means that, *on average*, cases with a
profile like this one were alive at five years 80% of the time in our reference
cohorts. It describes a group, not a certainty for any one patient, and it
carries statistical uncertainty (shown, where available, as a confidence
interval). The estimate is most reliable when the sample resembles the cases the
models were built on — the correct molecular subgroup, comparable assay quality,
and a similar treatment context.

---

> **Important.** These outputs are for **research and educational purposes only**.
> They are model-based, population-level estimates and are **not** a diagnosis, a
> clinical prognosis, or a basis for individual treatment decisions. They are
> intended to be interpreted by qualified researchers and clinicians alongside
> histopathology, established molecular classification, and full clinical context.
