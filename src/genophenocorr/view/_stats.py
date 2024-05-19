import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader
from genophenocorr.analysis import HpoMtcReport

from genophenocorr.model import Cohort


class StatsViewable:
    """
    Class to create a viewable object that is uses a Jinja2 template to create an HTML element
    for display in the Jupyter notebook.
    """

    def __init__(
            self,
            hpo_mtc_report: HpoMtcReport,
    ):
        """
        Args:
            hpo_mtc_report(HpoMtcReport): summary of heuristic term filtering procedure
        """
        self._hpo_mtc_filter_name = hpo_mtc_report.filter_method
        self._mtc_name = hpo_mtc_report.mtc_method
        self._results_map = hpo_mtc_report.skipped_terms_dict
        if not isinstance(hpo_mtc_report.skipped_terms_dict, dict):
            raise ValueError(f"hpo_mtc_report.skipped_terms_dict must be disctionary but was {type(hpo_mtc_report.skipped_terms_dict)}")
        self._term_count = hpo_mtc_report.n_terms_tested
        environment = Environment(loader=(PackageLoader('genophenocorr.view', 'templates')))
        self._cohort_template = environment.get_template("stats.html")

    def process(
            self,
            cohort: Cohort,
    ) -> str:
        """
        Create an HTML that should be shown with display(HTML(..)) of the ipython package.

        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook

        Returns:
            str: an HTML string with parameterized template for rendering
        """
        context = self._prepare_context(cohort)
        return self._cohort_template.render(context)

    def _prepare_context(
            self,
            cohort: Cohort,
    ) -> typing.Mapping[str, typing.Any]:
       
        skip_counts = list()
        sort_by_value = dict(sorted(self._results_map.items(), key=lambda item: item[1], reverse=True))
        for skipped, c in sort_by_value.items():
            skip_counts.append({"Skipped": skipped, "count": c})
        # The following dictionary is used by the Jinja2 HTML template
        return {
            "mtc_name": self._mtc_name,
            "hpo_mtc_filter_name": self._hpo_mtc_filter_name ,
            "total_hpo_count": self._term_count,
            "skipped_counts": skip_counts,
        }
