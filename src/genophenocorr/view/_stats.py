import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader

from genophenocorr.model import Cohort


class StatsViewable:
    """
    Class to create a viewable object that is uses a Jinja2 template to create an HTML element
    for display in the Jupyter notebook.
    """

    def __init__(
            self,
            filter_method_name: str,
            mtc_name:str,
            filter_results_map: typing.Dict[str,int],
            term_count: int ,
    ):
        """
        Args:
            filter_name(str): name of the method used to limit number of terms tested
            mtc_name(str): name of multiple testing correction procedure
            filter_results_map(Dict): Dictionary with counts and reasons for removed terms
            total_terms_tested(int): The total number of HPO terms tested after MTC limitation procedures
        """
        self._hpo_mtc_filter_name = filter_method_name
        self._mtc_name = mtc_name
        self._results_map = filter_results_map
        self._term_count = term_count
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
