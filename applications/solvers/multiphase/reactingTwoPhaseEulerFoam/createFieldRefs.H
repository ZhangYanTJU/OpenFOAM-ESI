phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();

const volScalarField& alpha1 = phase1;
const volScalarField& alpha2 = phase2;

volVectorField& U1 = phase1.URef();
surfaceScalarField& phi1 = phase1.phiRef();
const tmp<surfaceScalarField> talphaPhi1 = phase1.alphaPhi();
const auto& alphaPhi1 = talphaPhi1();

volVectorField& U2 = phase2.URef();
surfaceScalarField& phi2 = phase2.phiRef();
const tmp<surfaceScalarField> talphaPhi2 = phase2.alphaPhi();
const auto& alphaPhi2 = talphaPhi2();

surfaceScalarField& phi = fluid.phi();

rhoThermo& thermo1 = phase1.thermoRef();
rhoThermo& thermo2 = phase2.thermoRef();

volScalarField& rho1 = thermo1.rho();
const volScalarField& psi1 = thermo1.psi();

volScalarField& rho2 = thermo2.rho();
const volScalarField& psi2 = thermo2.psi();

const IOMRFZoneList& MRF = fluid.MRF();
fv::options& fvOptions = fluid.fvOptions();
