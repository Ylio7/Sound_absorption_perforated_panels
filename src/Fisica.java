
public class Fisica {

    final private double ViscoCinem = 1.85 * Math.pow(10, -5);
    final private double CteGas = 287.05;
    final private double ConduTerm = .0241;
    final private double CapCalorEspe = 1.01;
    final private double DensEspe = 1.293;
    final private double RelCalEspe = 1.402;
    final private double Pres1Atm = 101325;
    final private double resisflujo = 9100;
    final private double TempAire = 20;
    final private double PresionAtm = 1;
    final private double DensAire = (PresionAtm * Pres1Atm) / (CteGas * (273.15 + TempAire));
    final private double VelSonido = Math.sqrt(RelCalEspe * Pres1Atm / DensEspe) * Math.sqrt(1 + (TempAire / 273.15));
    final private double ZAire = VelSonido * DensAire;
    final private double TermInter1 = 2 * Math.PI / VelSonido;
    final private double TermInter2 = VelSonido / (2 * Math.PI);          //to be used...
    final private double RatioDensViscoAire = DensAire / ViscoCinem;

    // Constructor
    public Fisica() {
    }

    // Getters
    public double getViscoCinem() {
        return ViscoCinem;
    }

    public double getCteGas() {
        return CteGas;
    }

    public double getConduTerm() {
        return ConduTerm;
    }

    public double getCapCalorEspe() {
        return CapCalorEspe;
    }

    public double getDensEspe() {
        return DensEspe;
    }

    public double getRelCalEspe() {
        return RelCalEspe;
    }

    public double getPres1Atm() {
        return Pres1Atm;
    }

    public double getResisflujo() {
        return resisflujo;
    }

    public double getTempAire() {
        return TempAire;
    }

    public double getPresionAtm() {
        return PresionAtm;
    }

    public double getDensAire() {
        return DensAire;
    }

    public double getVelSonido() {
        return VelSonido;
    }

    public double getZAire() {
        return ZAire;
    }

    public double getTermInter1() {
        return TermInter1;
    }

    public double getTermInter2() {
        return TermInter2;
    }

    public double getRatioDensViscoAire() {
        return RatioDensViscoAire;
    }
}
