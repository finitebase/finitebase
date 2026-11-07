// Simcenter STAR-CCM+ macro: setting_up.java
// Written by Simcenter STAR-CCM+ 20.06.010
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.material.*;
import star.base.report.*;
import star.coupledflow.*;
import star.turbulence.*;
import star.vis.*;
import star.flow.*;
import star.kwturb.*;
import star.energy.*;
import star.walldistance.*;
import star.automation.*;









public class setting_up extends StarMacro {

  public void execute() {
    execute0();
    // E:/temp/blue/difusor/optimizer test/passage_simulation.sim
    execute1();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    ImportManager importManager_0 = 
      simulation_0.getImportManager();

    importManager_0.importMeshFiles(new StringVector(new String[] {resolvePath("E:\\temp\\blue\\difusor\\stator_passage.msh")}), NeoProperty.fromString("{\'FileOptions\': [{\'Sequence\': 44}]}"));

    FvRepresentation fvRepresentation_0 = 
      ((FvRepresentation) simulation_0.getRepresentationManager().getObject("Volume Mesh"));

    Region region_0 = 
      simulation_0.getRegionManager().getRegion("fluid");

    fvRepresentation_0.generateCompactMeshReport(new ArrayList<>(Arrays.<Region>asList(region_0)));

    simulation_0.getSceneManager().createGeometryScene("Geometry Scene", "Outline", "Surface", 1, null);

    Scene scene_0 = 
      simulation_0.getSceneManager().getScene("Geometry Scene 1");

    scene_0.resetCamera();

    CurrentView currentView_0 = 
      scene_0.getCurrentView();

    currentView_0.setInput(new DoubleVector(new double[] {0.01970784323899195, -0.022581376263033894, -8.001020119022817E-4}), new DoubleVector(new double[] {0.01970784323899195, -0.022581376263033894, 0.12458415228466418}), new DoubleVector(new double[] {0.0, 1.0, 0.0}), 0.04558541573243843, 0, 30.0);

    PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));

    physicsContinuum_0.enable(SingleComponentGasModel.class);

    physicsContinuum_0.enable(CoupledFlowModel.class);

    physicsContinuum_0.enable(IdealGasModel.class);

    physicsContinuum_0.enable(CoupledEnergyModel.class);

    physicsContinuum_0.enable(TurbulentModel.class);

    physicsContinuum_0.enable(RansTurbulenceModel.class);

    physicsContinuum_0.enable(KOmegaTurbulence.class);

    physicsContinuum_0.enable(SstKwTurbModel.class);

    physicsContinuum_0.enable(KwAllYplusWallTreatment.class);

    physicsContinuum_0.enable(SteadyModel.class);

    physicsContinuum_0.enable(GammaTransitionModel.class);

    Boundary boundary_0 = 
      region_0.getBoundaryManager().getBoundary("inlet");

    StagnationBoundary stagnationBoundary_0 = 
      ((StagnationBoundary) simulation_0.get(ConditionTypeManager.class).get(StagnationBoundary.class));

    boundary_0.setBoundaryType(stagnationBoundary_0);

    boundary_0.getConditions().get(FlowDirectionOption.class).setSelected(FlowDirectionOption.Type.COMPONENTS);

    FlowDirectionProfile flowDirectionProfile_0 = 
      boundary_0.getValues().get(FlowDirectionProfile.class);

    Units units_0 = 
      ((Units) simulation_0.getUnitsManager().getObject(""));

    flowDirectionProfile_0.getMethod(ConstantVectorProfileMethod.class).getQuantity().setComponentsAndUnits(0.951, -0.3091, 0.0, units_0);

    SupersonicPressureProfile supersonicPressureProfile_0 = 
      boundary_0.getValues().get(SupersonicPressureProfile.class);

    Units units_1 = 
      ((Units) simulation_0.getUnitsManager().getObject("Pa"));

    supersonicPressureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(294420.0, units_1);

    TotalPressureProfile totalPressureProfile_0 = 
      boundary_0.getValues().get(TotalPressureProfile.class);

    totalPressureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(360510.0, units_1);

    TotalTemperatureProfile totalTemperatureProfile_0 = 
      boundary_0.getValues().get(TotalTemperatureProfile.class);

    Units units_2 = 
      ((Units) simulation_0.getUnitsManager().getObject("K"));

    totalTemperatureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(443.1634, units_2);

    Boundary boundary_1 = 
      region_0.getBoundaryManager().getBoundary("outlet");

    MassFlowBoundary massFlowBoundary_0 = 
      ((MassFlowBoundary) simulation_0.get(ConditionTypeManager.class).get(MassFlowBoundary.class));

    boundary_1.setBoundaryType(massFlowBoundary_0);

    boundary_1.getConditions().get(MassFlowOption.class).setSelected(MassFlowOption.Type.MASS_FLUX);

    MassFluxProfile massFluxProfile_0 = 
      boundary_1.getValues().get(MassFluxProfile.class);

    Units units_3 = 
      ((Units) simulation_0.getUnitsManager().getObject("kg/m^2-s"));

    massFluxProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(-154.0, units_3);

    SupersonicPressureProfile supersonicPressureProfile_1 = 
      boundary_1.getValues().get(SupersonicPressureProfile.class);

    supersonicPressureProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(331840.0, units_1);

    Boundary boundary_2 = 
      region_0.getBoundaryManager().getBoundary("inlet_period_1");

    Boundary boundary_3 = 
      region_0.getBoundaryManager().getBoundary("inlet_period_2");

    BoundaryInterface boundaryInterface_0 = 
      simulation_0.getInterfaceManager().createBoundaryInterface(boundary_2, boundary_3, "Interface");

    Boundary boundary_4 = 
      region_0.getBoundaryManager().getBoundary("outlet_period_1");

    Boundary boundary_5 = 
      region_0.getBoundaryManager().getBoundary("outlet_period_2");

    BoundaryInterface boundaryInterface_1 = 
      simulation_0.getInterfaceManager().createBoundaryInterface(boundary_4, boundary_5, "Interface");

    boundaryInterface_0.getTopology().setSelected(InterfaceConfigurationOption.Type.PERIODIC);

    boundaryInterface_1.getTopology().setSelected(InterfaceConfigurationOption.Type.PERIODIC);

    Units units_4 = 
      ((Units) simulation_0.getUnitsManager().getObject("m"));

    physicsContinuum_0.getReferenceValues().get(MinimumAllowableWallDistance.class).setValueAndUnits(1.0E-7, units_4);

    physicsContinuum_0.getReferenceValues().get(ReferencePressure.class).setValueAndUnits(0.0, units_1);

    InitialPressureProfile initialPressureProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(InitialPressureProfile.class);

    initialPressureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(300000.0, units_1);

    StaticTemperatureProfile staticTemperatureProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(StaticTemperatureProfile.class);

    staticTemperatureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(400.0, units_2);

    VelocityProfile velocityProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(VelocityProfile.class);

    Units units_5 = 
      ((Units) simulation_0.getUnitsManager().getObject("m/s"));

    velocityProfile_0.getMethod(ConstantVectorProfileMethod.class).getQuantity().setComponentsAndUnits(0.0, -100.0, 0.0, units_5);












   CoupledImplicitSolver coupledImplicitSolver_0 = 
      ((CoupledImplicitSolver) simulation_0.getSolverManager().getSolver(CoupledImplicitSolver.class));

    CourantNumberController courantNumberController_0 = 
      coupledImplicitSolver_0.getCourantNumberController();

    AutomaticCourantNumberControl automaticCourantNumberControl_0 = 
      ((AutomaticCourantNumberControl) courantNumberController_0.getCourantNumberControl());

    automaticCourantNumberControl_0.setCFLMaximumValue(1000.0);

    automaticCourantNumberControl_0.setTargetAmgCycles(2);

    ExplicitRelaxationController explicitRelaxationController_0 = 
      coupledImplicitSolver_0.getExplicitRelaxationController();

    ConstantExplicitRelaxationControl constantExplicitRelaxationControl_0 = 
      ((ConstantExplicitRelaxationControl) explicitRelaxationController_0.getExplicitRelaxationControl());

    Units units_6 = 
      ((Units) simulation_0.getUnitsManager().getObject(""));

    constantExplicitRelaxationControl_0.getUrfQuantity().setValueAndUnits(0.1, units_6);















    currentView_0.setInput(new DoubleVector(new double[] {0.026802238651829593, -0.021803355098321042, -0.005035827592191761}), new DoubleVector(new double[] {0.026802238651829593, -0.021803355098321042, 0.15074682425647729}), new DoubleVector(new double[] {0.0, 1.0, 0.0}), 0.04558541573243843, 0, 30.0);

    StepStoppingCriterion stepStoppingCriterion_0 = 
      ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));

    
    IntegerValue integerValue_0 = 
      stepStoppingCriterion_0.getMaximumNumberStepsObject();

    integerValue_0.getQuantity().setValue(2000.0);















    ResidualPlot residualPlot_0 = 
      ((ResidualPlot) simulation_0.getPlotManager().getPlot("Residuals"));

    residualPlot_0.openInteractive();



    MassFlowAverageReport massFlowAverageReport_0 = 
      simulation_0.getReportManager().create("star.flow.MassFlowAverageReport");

    PrimitiveFieldFunction primitiveFieldFunction_0 = 
      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("AbsoluteTotalPressure"));

    massFlowAverageReport_0.setFieldFunction(primitiveFieldFunction_0);

    massFlowAverageReport_0.getParts().setQuery(null);

    massFlowAverageReport_0.getParts().setObjects(boundary_1);

    massFlowAverageReport_0.setSignedMassFlow(true);

    massFlowAverageReport_0.setSmooth(true);

    simulation_0.getMonitorManager().createMonitors(new ArrayList<>(Arrays.<Report>asList(massFlowAverageReport_0)), new PlotCreationInfo(PlotCreationInfo.CreatePlotChoice.SINGLE_PLOT, " Plot"));

    Cartesian2DPlot cartesian2DPlot_0 = 
      ((Cartesian2DPlot) simulation_0.getPlotManager().getPlot("Mass Flow Averaged 1 Monitor Plot"));

    cartesian2DPlot_0.openInteractive();

    PartDisplayer partDisplayer_0 = 
      ((PartDisplayer) scene_0.getDisplayerManager().getObject("Surface 1"));

    partDisplayer_0.setSurface(false);

    partDisplayer_0.setMesh(true);

    simulation_0.getSceneManager().createScalarScene("Scalar Scene", "Outline", "Scalar", null);

    Scene scene_1 = 
      simulation_0.getSceneManager().getScene("Scalar Scene 1");

    scene_1.resetCamera();

    scene_1.setPresentationName("pressure");

    ScalarDisplayer scalarDisplayer_0 = 
      ((ScalarDisplayer) scene_1.getDisplayerManager().getObject("Scalar 1"));

    PrimitiveFieldFunction primitiveFieldFunction_1 = 
      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("AbsolutePressure"));

    scalarDisplayer_0.getScalarDisplayQuantity().setFieldFunction(primitiveFieldFunction_1);

    scalarDisplayer_0.setUseIdealizations(true);

    simulation_0.getSceneManager().createScalarScene("Scalar Scene", "Outline", "Scalar", null);

    Scene scene_2 = 
      simulation_0.getSceneManager().getScene("Scalar Scene 1");

    scene_2.resetCamera();

    scene_2.setPresentationName("velocity");

    ScalarDisplayer scalarDisplayer_1 = 
      ((ScalarDisplayer) scene_2.getDisplayerManager().getObject("Scalar 1"));

    scalarDisplayer_1.setUseIdealizations(true);

    PrimitiveFieldFunction primitiveFieldFunction_2 = 
      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("Velocity"));

    VectorMagnitudeFieldFunction vectorMagnitudeFieldFunction_0 = 
      ((VectorMagnitudeFieldFunction) primitiveFieldFunction_2.getMagnitudeFunction());

    scalarDisplayer_1.getScalarDisplayQuantity().setFieldFunction(vectorMagnitudeFieldFunction_0);

    simulation_0.getSceneManager().createScalarScene("Scalar Scene", "Outline", "Scalar", null);

    Scene scene_3 = 
      simulation_0.getSceneManager().getScene("Scalar Scene 1");

    scene_3.resetCamera();

    scene_3.setPresentationName("k");

    ScalarDisplayer scalarDisplayer_2 = 
      ((ScalarDisplayer) scene_3.getDisplayerManager().getObject("Scalar 1"));

    scalarDisplayer_2.setUseIdealizations(true);

    PrimitiveFieldFunction primitiveFieldFunction_3 = 
      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("TurbulentKineticEnergy"));

    scalarDisplayer_2.getScalarDisplayQuantity().setFieldFunction(primitiveFieldFunction_3);

















   SimDriverWorkflow simDriverWorkflow_0 = 
      simulation_0.get(SimDriverWorkflowManager.class).createSimDriverWorkflow("Simulation Operations");

    ClearSolutionAutomationBlock clearSolutionAutomationBlock_0 = 
      (ClearSolutionAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.automation.ClearSolutionAutomationBlock", "Clear Solution");

    InitializeSolutionAutomationBlock initializeSolutionAutomationBlock_0 = 
      (InitializeSolutionAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.automation.InitializeSolutionAutomationBlock", "Initialize Solution");

    SolvePhysics solvePhysics_0 = 
      (SolvePhysics) simDriverWorkflow_0.getBlocks().createBlock("star.common.SolvePhysics", "Solve Physics");

    solvePhysics_0.getSimulationObjects().setQuery(null);

    PhysicsContinuum physicsContinuum_1 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));

    solvePhysics_0.getSimulationObjects().setObjects(physicsContinuum_1);

    ExportAutomationBlock exportAutomationBlock_0 = 
      (ExportAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.vis.ExportAutomationBlock", "Export");

    exportAutomationBlock_0.getSources().setQuery(null);

    //Cartesian2DPlot cartesian2DPlot_0 = 
     // ((Cartesian2DPlot) simulation_0.getPlotManager().getPlot("Mass Flow Averaged 1 Monitor Plot"));

    exportAutomationBlock_0.getSources().setObjects(cartesian2DPlot_0);

    PlotExportSettings plotExportSettings_0 = 
      exportAutomationBlock_0.getPlotSettings();

    plotExportSettings_0.setExportType(PlotExportType.CSV);



    simulation_0.get(SimDriverWorkflowManager.class).setSelectedWorkflow(simDriverWorkflow_0);





    simDriverWorkflow_0.execute();


    System.exit(0);



















   









  }

  private void execute1() {
  }
}
