import hpotk
import pandas as pd

from gpsea.analysis.pcats import HpoTermAnalysisResult


def summarize_hpo_analysis(
    hpo: hpotk.MinimalOntology,
    result: HpoTermAnalysisResult,
) -> pd.DataFrame:
    """
    Create a dataframe with counts, frequencies, and p values for the tested HPO terms.

    The HPO terms that were not tested will *not* be included in the frame.

    :param hpo: HPO data.
    :param result: the HPO term analysis results to show.
    """
    assert isinstance(result, HpoTermAnalysisResult)

    # Row index: a list of tested HPO terms
    pheno_idx = pd.Index(result.phenotypes)
    # Column index: multiindex of counts and percentages for all genotype predicate groups
    gt_idx = pd.MultiIndex.from_product(
        iterables=(result.gt_predicate.get_categories(), ("Count", "Percent")),
        names=(result.gt_predicate.get_question_base(), None),
    )

    # We'll fill this frame with data
    df = pd.DataFrame(index=pheno_idx, columns=gt_idx)

    for ph_predicate, count in zip(result.pheno_predicates, result.all_counts):
        # Sum across the phenotype categories (collapse the rows).
        gt_totals = count.sum()

        for gt_cat in count.columns:
            cnt = count.loc[ph_predicate.present_phenotype_category, gt_cat]
            total = gt_totals[gt_cat]
            df.loc[ph_predicate.phenotype, (gt_cat, "Count")] = f"{cnt}/{total}"
            pct = 0 if total == 0 else round(cnt * 100 / total)
            df.loc[ph_predicate.phenotype, (gt_cat, "Percent")] = f"{pct}%"

    # Add columns with p values and corrected p values (if present)
    p_val_col_name = "p values"
    corrected_p_val_col_name = "Corrected p values"
    if result.corrected_pvals is not None:
        df.insert(df.shape[1], ("", corrected_p_val_col_name), result.corrected_pvals)
    df.insert(df.shape[1], ("", p_val_col_name), result.pvals)

    # Format the index values: `HP:0001250` -> `Seizure [HP:0001250]` if the index members are HPO terms
    # or just use the term ID CURIE otherwise (e.g. `OMIM:123000`).
    labeled_idx = df.index.map(lambda term_id: format_term_id(hpo, term_id))

    # Last, sort by corrected p value or just p value
    df = df.set_index(labeled_idx)
    # and only report the tested HPO terms
    with_p_value = df[("", p_val_col_name)].notna()
    if result.corrected_pvals is not None:
        return df.sort_values(by=[("", corrected_p_val_col_name), ("", p_val_col_name)]).loc[with_p_value]
    else:
        return df.sort_values(by=("", p_val_col_name)).loc[with_p_value]


def format_term_id(
    hpo: hpotk.MinimalOntology,
    term_id: hpotk.TermId,
) -> str:
    """
    Format a `term_id` as a `str`. HPO term ID is formatted as `<name> [<term_id>]` whereas other term IDs
    are formatted as CURIEs (e.g. `OMIM:123000`).
    """
    if term_id.prefix == "HP":
        term_name = hpo.get_term_name(term_id)
        return f"{term_name} [{term_id.value}]"
    else:
        return term_id.value
