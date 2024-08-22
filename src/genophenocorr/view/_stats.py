import typing

from jinja2 import Environment, PackageLoader
from genophenocorr.analysis import HpoMtcReport


class MtcStatsViewer:
    """
    `MtcStatsViewer` uses a Jinja2 template to create an HTML element for showing in the Jupyter notebook
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
        Create an HTML that should be shown with `display(HTML(..))` of the IPython package.

        Args:
            hpo_mtc_report (HpoMtcReport): summary of heuristic term filtering procedure

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

        n_skipped = 0
        reason_to_count = list()
        for reason, count in sorted(results_map.items(), key=lambda item: item[1], reverse=True):
            reason_to_count.append({"reason": reason, "count": count})
            n_skipped += count

        n_tested = hpo_mtc_report.n_terms_before_filtering - n_skipped

        # The following dictionary is used by the Jinja2 HTML template
        return {
            "mtc_name": hpo_mtc_report.mtc_method,
            "hpo_mtc_filter_name": hpo_mtc_report.filter_method,
            "skipped_hpo_count": n_skipped,
            "tested_hpo_count": n_tested,
            "total_hpo_count": hpo_mtc_report.n_terms_before_filtering,
            "reason_to_count": reason_to_count,
        }
