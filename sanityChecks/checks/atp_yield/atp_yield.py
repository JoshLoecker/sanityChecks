import re
import pandas as pd
from enum import Enum
from cobra import Model, Reaction, Metabolite, Solution

from args import Arguments


class _Metabolites(Enum):
    OXYGEN = "EX_o2[e]"
    WATER = "EX_h2o[e]"
    CARBON_DIOXIDE = "EX_co2[e]"
    GLUCOSE = "EX_glc_D[e]"
    GLUTAMINE = "EX_gln_L[e]"
    FRUCTOSE = "EX_fru[e]"
    BUTYRATE = "EX_but[e]"
    CAPROIC_ACID = "EX_caproic[e]"
    OCTANOATE = "EX_octa[e]"
    DECANOATE = "EX_dca[e]"
    LAUREATE = "EX_ddca[e]"
    TETRADECANOATE = "EX_ttdca[e]"
    HEXADECANOATE = "EX_hdca[e]"
    OCTADECANOATE = "EX_ocdca[e]"
    ARACHIDATE = "EX_arach[e]"
    BEHENIC_ACID = "EX_docosac[e]"
    LIGNOCERATE = "EX_lgnc[e]"
    CEROTATE = "EX_hexc[e]"


class _MetaboliteName(Enum):
    EX_o2 = "Extracellular oxygen exchange"
    EX_h2o = "Extracellular water exchange"
    EX_co2 = "Extracellular carbon dioxide exchange"
    EX_glc_D = "Extracellular glucose exchange"
    EX_gln_L = "Extracellular glutamine exchange"
    EX_fru = "Extracellular fructose exchange"
    EX_but = "Extracellular butyrate exchange"
    EX_caproic = "Extracellular caproic exchange"
    EX_octa = "Extracellular octanoate exchange"
    EX_dca = "Extracellular decanoate exchange"
    EX_ddca = "Extracellular laureate exchange"
    EX_ttdca = "Extracellular tetradecanoate exchange"
    EX_hdca = "Extracellular hexadecanoate exchange"
    EX_ocdca = "Extracellular octadecanoate exchange"
    EX_arach = "Extracellular arachidate exchange"
    EX_docosac = "Extracellular behenic exchange"
    EX_lgnc = "Extracellular lignocerate exchange"
    EX_hexc = "Extracellular cerotate exchange"


class ATPYield:
    def __init__(self, arguments: Arguments):
        self.model: Model = arguments.model
        self.solver: str = arguments.solver
        self.is_human: bool = arguments.is_human
        self._atp_synthase: str = "ATPS4mi" if self.is_human else "ATPS4m"
        self.original_closed_model: Model
        self.min_flux: float = 1e-6  # Any flux lower than this will get set to 0
        
        self.results_table: pd.DataFrame = pd.DataFrame(
            data={},
            index=[
                f"{self.model.name}: ATP yield",
                f"{self.model.name}: ATPS4m yield",
                "Theoretical"
            ]
        )
        
        self.__post_init__()
    
    def __post_init__(self):
        self.model.solver = self.solver
        if self.model.name is None:
            self.model.name = "Model"
        
        # Add ATPS4mi reaction if model is human, otherwise add ATPS4m
        self.model.add_reactions([
            Reaction(
                id=self._atp_synthase,
                name="ATP synthase (mitochondrial)",
                lower_bound=-1000,
                upper_bound=1000,
            )
        ])
        
        self._harmonize_reactions()
        self._close_model()
        self.test_atp_yield()
    
    def _harmonize_reactions(self):
        """
        We want to make sure all reactions have the same usage of brackets instead of parentheses.
        
        EXAMPLE:
        Initial: 10fthf6glu(m) --> 10fthf6glu(c)
        Result: 10fthf6glu[m] --> 10fthf6glu[c]
        :return:
        """
        reaction: Reaction
        for reaction in self.model.reactions:
            reaction.reaction = re.sub(pattern=r"\(", repl="[", string=reaction.reaction)
            reaction.reaction = re.sub(pattern=r"\(", repl="]", string=reaction.reaction)
        
        metabolite: Metabolite
        for metabolite in self.model.metabolites:
            metabolite.id = re.sub(pattern=r"\(", repl="[", string=metabolite.id)
            metabolite.id = re.sub(pattern=r"\)", repl="]", string=metabolite.id)
    
    def _close_model(self) -> None:
        """
        This function will set the lower bound of all exchange and sink reactions to ensure that only those metabolites that are supposed to be taken up are indeed taken up.
        """
        closed_model = self.model.copy()
        
        # Find all reactions based on their abbreviation
        exchange_reactions: list[Reaction] = [
            reaction for reaction in closed_model.reactions if
            str(reaction.id).upper().startswith("EX_")
        ]
        demand_reactions: list[Reaction] = [
            reaction for reaction in closed_model.reactions if
            str(reaction.id).upper().startswith("DM_")
        ]
        sink_reactions: list[Reaction] = [
            reaction for reaction in closed_model.reactions if
            str(reaction.id).upper().startswith("SINK_") or str(reaction.id).upper().startswith("SK_")
        ]
        
        # Get all reactions that have "bioma" in their name
        biomass_reactions: list[Reaction] = [
            reaction for reaction in closed_model.reactions if
            "bioma" in str(reaction.id).lower()
        ]
        
        # Select all all reactions that contain only one non-zero entry
        selected_exchange = [
            reaction for reaction in closed_model.reactions if
            sum(reaction.metabolites.values()) == 1
        ]
        
        # Collect exchange, demand, and sinks using Cobra's function (in case some are different)
        cobra_exchanges: list[Reaction] = self.model.exchanges
        cobra_demands: list[Reaction] = self.model.demands
        cobra_sinks: list[Reaction] = self.model.sinks
        
        # Combine the above lists into one list containing unique reactions
        model_reactions: set[Reaction] = {
            *exchange_reactions,
            *demand_reactions,
            *sink_reactions,
            *biomass_reactions,
            *selected_exchange,
            *cobra_exchanges,
            *cobra_demands,
            *cobra_sinks
        }
        
        # For each reaction present in `model_reaction`, find that reaction `closed_model`, and set its lower bound to 0
        reaction: Reaction
        for reaction in model_reactions:
            closed_model_reaction: Reaction = closed_model.reactions.get_by_id(reaction.id)
            closed_model_reaction.lower_bound = 0
        
        # Now, set all reactions in `selected_exchange` to have an upper bound of 1000 (infinity)
        for reaction in selected_exchange:
            closed_model_reaction = closed_model.reactions.get_by_id(reaction.id)
            closed_model_reaction.upper_bound = 1000
        
        # Add `DM_apt_c` to the model and set it to be the objective
        closed_model.add_metabolites([
            Metabolite(
                id="atp[c]",
                name="Cystolic APT",
                compartment="c"
            )
        ])
        closed_model.add_boundary(closed_model.metabolites.get_by_id("atp[c]"), type="demand")
        closed_model.objective = closed_model.reactions.get_by_id("DM_atp[c]")
        self.original_closed_model = closed_model.copy()
    
    def _set_carbon_source(
        self,
        source: _Metabolites,
        allow_oxygen: bool,
        column_name: str = None,
    ) -> Model:
        """
        This function will set up a model with the following values:
            - Oxygen
                - If `allow_oxygen` is True, the lower bound will be set to -1000
                - If `allow_oxygen` if False, the lower bound AND upper bound will be set to 0
            - Water
                - The upper will always be set to 1000
                - The lower bound will always be set to -1000
            - Carbon Dioxide: The upper bound will always be set to 1000
            - `source`: The upper AND lower bound will always be set to -1
        :param source: An item from the CarbonSource enum
        :param column_name: The name of the column that will be added to the results. If none, defaults to "`metabolite` - (an)aeorobic"
        :param allow_oxygen: A boolean value indicating whether to allow oxygen to be used in the model
        :return: A model with the above indicated items set
        """
        if column_name is None:
            metabolite = source.value.split("_")[1]  # example: extract "glc" from "EX_glc_D[e]"
            column_name = f"{metabolite} - aerobic" if allow_oxygen else f"{metabolite} - anaerobic"
        
        model_copy: Model = self.original_closed_model.copy()
        model_copy.add_metabolites([
            Metabolite(id="o2[e]", name="Extracellular oxygen", compartment="e"),
            Metabolite(id="h2o[e]", name="Extracellular water", compartment="e"),
            Metabolite(id="co2[e]", name="Extracellular carbon dioxide", compartment="e")
        ])
        
        model_copy.add_reactions([
            Reaction(
                id=_Metabolites.OXYGEN.value,
                name="Extracellular oxygen exchange",
                lower_bound=-1000 if allow_oxygen else 0,
                upper_bound=1000 if allow_oxygen else 0
            ),
            Reaction(
                id=_Metabolites.WATER.value,
                name="Extracellular water exchange",
                lower_bound=-1000,
                upper_bound=1000
            ),
            Reaction(
                id=_Metabolites.CARBON_DIOXIDE.value,
                name="Extracellular carbon dioxide exchange",
                lower_bound=-1000,
                upper_bound=1000
            )
        ])
        
        # Test if the current carbon source is already in the model
        # If it is, set its lower and upper bound to -1
        try:
            model_copy.reactions.get_by_id(source.value).lower_bound = -1
            model_copy.reactions.get_by_id(source.value).upper_bound = -1
        except KeyError:
            reaction_name = source.value.split("[")[0]  # Get "EX_o2" from "EX_o2[e]"
            model_copy.add_reactions([
                Reaction(
                    id=source.value,
                    name=_MetaboliteName[reaction_name].value,
                    lower_bound=-1,
                    upper_bound=-1
                )
            ])
        
        # Compute FBA
        model_copy.objective = model_copy.reactions.get_by_id("DM_atp[c]")
        FBA: Solution = model_copy.optimize("maximize")
        
        # Determine which ATP synthase to use (ATPs4m or ATPs4mi)
        print("\n\nHERE")
        print(FBA.fluxes.index.tolist())
        
        if len(FBA.fluxes) > 0:
            # For each flux value in FBA.flux, if it is less than `self.min_flux`, set it to 0
            FBA.fluxes.apply(
                lambda flux_value: 0 if flux_value < self.min_flux else flux_value
            )
            self.results_table.insert(
                loc=len(self.results_table.columns),
                column=column_name,
                value=[
                    FBA.objective_value,  # Row 1
                    FBA.fluxes[self._atp_synthase],  # Row 2, ATPS4m
                    "31",  # Row 3, theoretical yield
                ]
            )
        
        print(self.results_table)
        exit(1)
        
        return model_copy
    
    def test_atp_yield(self):
        """
        This function will the ATP yield of the model using the following carbon sources/conditions:
            - Glucose:          [Unlimited oxygen, No oxygen]
            - Glutamine:        [Unlimited oxygen, No oxygen]
            - Fructose:         [Unlimited oxygen, No oxygen]
            - Butyrate:         [Unlimited oxygen, No oxygen]
            - Caproic acid:     [Unlimited oxygen, No oxygen]
            - Octanoate:        [Unlimited oxygen, No oxygen]
            - Decanoate:        [Unlimited oxygen, No oxygen]
            - Laureate:         [Unlimited oxygen, No oxygen]
            - Tetradecanoate:   [Unlimited oxygen, No oxygen]
            - Hexadecanoate:    [Unlimited oxygen, No oxygen]
            - Octadecanoate:    [Unlimited oxygen, No oxygen]
            - Arachidate:       [Unlimited oxygen, No oxygen]
            - Behenic:          [Unlimited oxygen, No oxygen]
            - Lignocerate:      [Unlimited oxygen, No oxygen]
            - Cerotate:         [Unlimited oxygen, No oxygen]
        :return:
        """
        # Test for APT yield with 1 mol/gdw/hr glucose and unlimited oxygen
        
        closed_model = self.original_closed_model.copy()
        
        # Glucose
        glucose_and_oxygen = self._set_carbon_source(source=_Metabolites.GLUCOSE, allow_oxygen=True)
        
        # glucose_no_oxygen = self._set_carbon_source(source=_Metabolites.GLUCOSE, allow_oxygen=False)
