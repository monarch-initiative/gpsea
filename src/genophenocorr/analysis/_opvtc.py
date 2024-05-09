
import typing
import pandas as pd
import hpotk
from statsmodels.stats import multitest
from collections import defaultdict, deque

from ._stats import run_fisher_exact
from .predicate import BooleanPredicate,PolyPredicate
from .predicate.phenotype import PhenotypePredicateFactory
from genophenocorr.model import Patient




class OntologyPValueTestChooser:
    """
        `OntologyPValueTestChooser` decides which P-Values tests to perform and which to skip

        The class uses a number of heuristics and rules of thumb.
        1. Do not perform a test if there are less than min_annots annotations (default 2). That is, do not perform a
            test if an HPO term is used only once for annotation
        2. Do not perform a test if the counts in the genotype categories do not even have nominal statistical power
            for i in range(2,6):
                for j in range(2,6):
                    # create a table with the most extreme possible counts. If this is not significant, then
                    # other counts also cannot be
                    two_by_two_table = [[i, 0], [0, j]]
                    oddsr, p = scipy.stats.fisher_exact(two_by_two_table, alternative='two-sided')
            This code reveals that (2,4), (4,2), (2,3), (3,3), (2,2) and (3,2) are not powered so we can always skip them
            ## TODO -- similar calculate for 3x2
    """
    def __init__(self, pheno_predicate_factory, patients: typing.Iterable[Patient],
                phenotypic_features: typing.Iterable[hpotk.TermId],
                predicate: PolyPredicate,
                hpo: hpotk.MinimalOntology,
                 second_level_pval_threshold=1e-5,
                 min_annots:int=2):
        """



        """
        self._phenopredicate_factory = pheno_predicate_factory
        self._patients = patients
        self._phenotypic_features = phenotypic_features
        self._predicate = predicate
        self._hpo = hpo
        self._min_annots = min_annots
        self._SECONDLEVEL_PVAL_THRES = second_level_pval_threshold
        # get all HPO terms in the graph that is induced by the phenotypic observations
        self._hpo_term_id_set = set()
        self._hpo_leaf_term_id_set = set()
        for pf in self._phenotypic_features:
            self._hpo_term_id_set.add(pf.id)
            # get ancestors (data type: Iterator)
            ancs = self._hpo.graph.get_ancestors(pf.id, False)
            ancs = set(ancs)
            self._hpo_term_id_set.update(ancs)
        print(f"[INFO] We found a total of {len(self._hpo_term_id_set)} unique directly and indirectly annotated HPO terms")
        # get all of the leaf terms
        for pf in self._phenotypic_features:
            # get descendants not including original term
            descs = self._hpo.graph.get_descendants(pf.id, False)
            descs = set(descs)
            # if no descendant is in the set of annotated terms, then the term must be a leaf term
            if len(descs.intersection(self._hpo_term_id_set)) == 0:
                self._hpo_leaf_term_id_set.add(pd.id)
        print(f"[INFO] We found a total of {len(self._hpo_leaf_term_id_set)} HPO terms as leaves of the graph")
        # The following numbers of total observations in the genotype groups can never be signficant so we just skip them
        # see above explanation
        self._powerless_set = {(2, 4), (4, 2), (2, 3), (3, 3), (2, 2), (3, 2)}
        # thus if the total count is 6 or less, there is no power
        ## Make an ordering of the HPO ids such that we start with the leaves and work our way upward, without skipping a term
        self._ordered_term_list = list()
        queue = deque()
        for hpo_id in self._hpo_leaf_term_id_set:
            queue.append(hpo_id)
        while len(queue) > 0:
            hpo_id = queue.popleft()
            self._ordered_term_list.append(hpo_id)
            parents = self._hpo.graph.get_parents(hpo_id, False)
            for p in parents:
                queue.append(p)
        # now, self._ordered_term_list is ordered from leaves to root
        # Finally, collect sets of top-level and second-level terms
        self._PHENO_ROOT_ID = "HP:0000118"
        self._top_level_terms = set()  # e.g., Abnormality of the cardiovascular system (HP:0001626)
        self._second_level_terms = set()
        for id in hpo.graph.get_children(self._PHENO_ROOT_ID, False):
            self._top_level_terms.add(id)
        for top_level_id in self._top_level_terms:
            for id in hpo.graph.get_children(top_level_id, False):
                self._second_level_terms.add(id)




    def _do_2by2_pval_test(self, geno_predicate:BooleanPredicate):
        """
        This method decides which P-value tests to perform based on a series of heuristics that do not
        test terms if they are impossible or very unlikely to provide an interesting result.
        :param geno_predicate: genotype predictate, e.g., Is this a nonsense variant?
        :type geno_predicate: BooleanPredicate
        """
        tested_pf = defaultdict(tuple) # key is an HP id, value is a tuple with counts, i.e.,
        pfs_with_pvals = defaultdict(float)
        # (geno_yes_pheno_no,geno_yes_pheno_yes,geno_no_pheno_no,geno_no_pheno_yes)
        refused_to_test = dict() ## Probably replace strings with an ENUM
        for hpo_id in self._ordered_term_list:
            if not self._hpo.graph.is_descendant_of(hpo_id, self._PHENO_ROOT_ID):
                refused_to_test[hpo_id] = "Not a phenotype term"
                continue ## only test phenotype annotations
            if hpo_id in self._top_level_terms:
                refused_to_test[hpo_id] = "Do not test top-level (general) terms"
                continue  ## only test phenotype annotations
            pheno_predicate = self._pheno_predicate_factory.get_predicate(hpo_id)
            total_count = 0
            geno_yes_pheno_no = 0
            geno_yes_pheno_yes = 0
            geno_no_pheno_no = 0
            geno_no_pheno_yes = 0
            for patient in self._patients:
                pheno_cat = pheno_predicate.test(patient)
                geno_cat = geno_predicate.test(patient)
                if pheno_cat is not None and geno_cat is not None:
                    if pheno_cat == BooleanPredicate.YES:
                        if geno_cat == BooleanPredicate.YES:
                            geno_yes_pheno_yes +=1
                        elif geno_cat == BooleanPredicate.NO:
                            geno_no_pheno_yes += 1
                        else:
                            raise ValueError(f"Unrecognized genocat {geno_cat}")
                    elif pheno_cat == BooleanPredicate.NO:
                        if geno_cat == BooleanPredicate.YES:
                            geno_yes_pheno_no +=1
                        elif geno_cat == BooleanPredicate.NO:
                            geno_yes_pheno_no += 1
                        else:
                            raise ValueError(f"Unrecognized geno_cat {geno_cat}")
                    else:
                        raise ValueError(f"Unrecognized pheno_cat {pheno_cat}")
                    total_count += 1
            if total_count <= 6:
                refused_to_test[hpo_id] = "Counts too low to have power"
                continue
            if (geno_yes_pheno_yes + geno_no_pheno_yes) <= 1:
                refused_to_test[hpo_id] = "Do not test if only a single positive observation of an HPO term"
                continue
            # If we get here, then the term has more than on positive annotation in our cohort and the
            # are enough counts to potentially achieve statistical power. Now we check if the term has a descendent that
            # was already tested and that has the same counts
            counts_tuple = (geno_yes_pheno_no,geno_yes_pheno_yes,geno_no_pheno_no,geno_no_pheno_yes)
            children = self._hpo.graph.get_children(hpo_id, False)
            for c in children:
                if c in tested_pf:
                    if tested_pf[c] == counts_tuple:
                        refused_to_test[hpo_id] = "Child term with same counts previously tested"
                        continue
            ## Heuristic -- if a child of a second level term was very significant, then do not test the
            # second level term, because it is unlikely to add much insight
            # A first level term is like Abnormality of the digestive system (HP:0025031)
            # A second level term is like Abnormality of digestive system physiology (HP:0025032)
            # A specific term in this hierarchy would be like Esophageal diverticulum HP:0100628
            if hpo_id in self._second_level_terms:
                doTest = True
                desc = self._hpo.graph.get_descendants(hpo_id, False)
                for d in desc:
                    if d in pfs_with_pvals and pfs_with_pvals[d] < self._SECONDLEVEL_PVAL_THRES:
                        refused_to_test[hpo_id] = "Highly significant descendant of second-level term"
                        doTest = False
                        break
                if not doTest:
                    continue
            # if we get here, then this term id has passed all of our heuristics, so we go ahead and test
            p_val = run_fisher_exact([[geno_yes_pheno_yes, geno_yes_pheno_no], [geno_no_pheno_yes, geno_no_pheno_no]])
            pfs_with_pvals[hpo_id] = p_val



