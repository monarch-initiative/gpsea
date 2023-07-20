from ._cohort_data import Cohort
from hpotk import Ontology


class CohortViewer:
    def __init__(self, cohort, hpo):
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort argument must be a Cohort object but we got a {type(cohort)} object")
        if not isinstance(hpo, Ontology):
            raise ValueError(f"cohort argument must be an Ontology but we got a {type(cohort)} object")
        self._cohort = cohort
        self._hpo = hpo

    @staticmethod
    def html_row(items):
        line = "</td><td>".join(items)
        return f"<tr><td>{line}</td></tr>"

    @staticmethod
    def anchor(url, payload):
        return f"<a href=\"{url}\">{payload}</a>"

    @staticmethod
    def hpo_anchor(hpo_id, label):
        url = f"https://hpo.jax.org/app/browse/term/{hpo_id}"
        return CohortViewer.anchor(url, label)

    def hpo_term_counts_table(self, min_count=1):
        rows = list()
        rows.append("<table>")
        header_items = ["HPO Term", "Count"]
        rows.append(CohortViewer.html_row(header_items))
        hpo_list = self._cohort.list_all_phenotypes()
        min_terms = set()
        for hpo_tuple in hpo_list:
            row_items = []
            hpo_id = hpo_tuple[0]
            count = hpo_tuple[1]
            term = self._hpo.get_term(hpo_id)
            if term is None:
                hpo_label = "n/a"
            else:
                hpo_label = term.name
            if count >= min_count:
                row_items.append(CohortViewer.hpo_anchor(hpo_id, hpo_label))
                row_items.append(str(count))
                rows.append(CohortViewer.html_row(row_items))
            else:
                min_terms.add(CohortViewer.hpo_anchor(hpo_id, hpo_label))
        rows.append("</table>")
        if len(min_terms) > 0:
            term_str = "; ".join(min_terms)
            rows.append(f"<p>Additionally, the following terms were observed not more than {min_count} times:</p>")
            rows.append(f"<p>{term_str}.</p>")
        return "\n".join(rows)
