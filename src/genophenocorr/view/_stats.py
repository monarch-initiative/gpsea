import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader
from genophenocorr.analysis import HpoMtcReport

from genophenocorr.model import Cohort


class StatsViewable:
    """
    `StatsViewable` uses a Jinja2 template to create an HTML element for showing in the Jupyter notebook
    or for writing into a standalone HTML file.
    """

    def __init__(self):
        environment = Environment(loader=(PackageLoader('genophenocorr.view', 'templates')))
        self._cohort_template = environment.get_template("stats.html")

    def process(
            self,
            hpo_mtc_report: HpoMtcReport,
    ) -> str:
        """
        Create an HTML that should be shown with `display(HTML(..))` of the ipython package.

        Args:
            hpo_mtc_report(HpoMtcReport): summary of heuristic term filtering procedure

        Returns:
            str: an HTML string with parameterized template for rendering or writing into a standalone HTML file.
        """
        context = self._prepare_context(hpo_mtc_report)
        return self._cohort_template.render(context)

    @staticmethod
    def _prepare_context(
            hpo_mtc_report: HpoMtcReport,
    ) -> typing.Mapping[str, typing.Any]:
        results_map = hpo_mtc_report.skipped_terms_dict
        if not isinstance(results_map, dict):
            raise ValueError(
                "hpo_mtc_report.skipped_terms_dict must be dictionary "
                f"but was {type(hpo_mtc_report.skipped_terms_dict)}"
            )
        skip_counts = list()
        sort_by_value = dict(sorted(results_map.items(), key=lambda item: item[1], reverse=True))
        for skipped, c in sort_by_value.items():
            skip_counts.append({"Skipped": skipped, "count": c})
        # The following dictionary is used by the Jinja2 HTML template
        return {
            "mtc_name": hpo_mtc_report.mtc_method,
            "hpo_mtc_filter_name": hpo_mtc_report.filter_method,
            "total_hpo_count": hpo_mtc_report.n_terms_tested,
            "skipped_counts": skip_counts,
        }
