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
    # Column index: a summary of the counts and percentages for all genotype predicate groups
    cat_names = list(cat.name for cat in result.gt_clf.get_categories())
    gt_idx = pd.Index(
        data=cat_names,
        name=result.gt_clf.variable_name,
    )

    # We'll fill this frame with data
    df = pd.DataFrame(index=pheno_idx, columns=gt_idx)

    for ph_clf, count in zip(result.pheno_clfs, result.all_counts):
        # Sum across the phenotype categories (collapse the rows).
        gt_totals = count.sum()

        for gt_cat in count.columns:
            cnt = count.loc[ph_clf.present_phenotype_category, gt_cat]
            total = gt_totals[gt_cat]
            pct = 0 if total == 0 else round(cnt * 100 / total)
            df.loc[ph_clf.phenotype, gt_cat.name] = f"{cnt}/{total} ({pct}%)"

    # Add columns with p values and corrected p values (if present)
    p_val_col_name = "p values"
    corrected_p_val_col_name = "Corrected p values"
    if result.corrected_pvals is not None:
        df.insert(df.shape[1], corrected_p_val_col_name, result.corrected_pvals)
    df.insert(df.shape[1], p_val_col_name, result.pvals)

    # Format the index values: `HP:0001250` -> `Seizure [HP:0001250]` if the index members are HPO terms
    # or just use the term ID CURIE otherwise (e.g. `OMIM:123000`).
    labeled_idx = df.index.map(lambda term_id: format_term_id(hpo, term_id))

    # Last, sort by corrected p value or just p value
    df = df.set_index(labeled_idx)
    # and only report the tested HPO terms
    with_p_value = df[p_val_col_name].notna()

    sort_columns = []
    if result.corrected_pvals is not None:
        sort_columns.append(corrected_p_val_col_name)
        sort_columns.append(p_val_col_name)
        sort_columns.extend(cat_names)
    else:
        sort_columns.append(p_val_col_name)
        sort_columns.extend(cat_names)

    return df.sort_values(by=sort_columns).loc[with_p_value]


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
