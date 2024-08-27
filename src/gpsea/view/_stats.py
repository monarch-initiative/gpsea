import typing

from collections import Counter

from jinja2 import Environment, PackageLoader

from gpsea.analysis.pcats import HpoTermAnalysisResult


class MtcStatsViewer:
    """
    `MtcStatsViewer` uses a Jinja2 template to create an HTML element for showing in the Jupyter notebook
    or for writing into a standalone HTML file.
    """

    def __init__(self):
        environment = Environment(loader=(PackageLoader('gpsea.view', 'templates')))
        self._cohort_template = environment.get_template("stats.html")

    def process(
        self,
        result: HpoTermAnalysisResult,
    ) -> str:
        """
        Create an HTML to present MTC part of the :class:`HpoTermAnalysisResult`.

        Use the `display(HTML(..))` functions of the IPython package.

        Args:
            result (HpoTermAnalysisResult): the result to show

        Returns:
            str: an HTML string with parameterized template for rendering or writing into a standalone HTML file.
        """
        assert isinstance(result, HpoTermAnalysisResult)
        context = self._prepare_context(result)
        return self._cohort_template.render(context)

    @staticmethod
    def _prepare_context(
        report: HpoTermAnalysisResult,
    ) -> typing.Mapping[str, typing.Any]:
        counts = Counter()
        for result in report.mtc_filter_results:
            if result.is_filtered_out():
                counts[result.reason] += 1

        n_skipped = 0
        reason_to_count = list()
        for reason, count in sorted(counts.items(), key=lambda item: item[1], reverse=True):
            reason_to_count.append({"reason": reason, "count": count})
            n_skipped += count

        n_all = len(report.phenotypes)
        n_tested = n_all - n_skipped

        # The following dictionary is used by the Jinja2 HTML template
        return {
            "mtc_name": report.mtc_name,
            "hpo_mtc_filter_name": report.mtc_filter_name,
            "skipped_hpo_count": n_skipped,
            "tested_hpo_count": n_tested,
            "total_hpo_count": n_all,
            "reason_to_count": reason_to_count,
        }
