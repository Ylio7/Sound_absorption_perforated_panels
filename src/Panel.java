
public class Panel {

    final private double EspesorPanel;
    final private double DistanciaRepeticion;
    final private double RadioPerf;
    final private double cavidad;
    final private double espesorAire;
    final private double areaAbierta;


    final private double endCorr;
    final private double freq[];
    final private double w[];

    final private Complejo[] Kg;
    final private Complejo[] CotKaireTaire;
    final private Complejo CotKaTa[];
    final private Complejo J[];
    final private Complejo negJ[];
    final private Complejo Zc[];
    final private Complejo ZAire_c[];

    final private Fisica F = new Fisica();

    public Panel(double EspesorPanel, double DistanciaRepeticion, double RadioPerf, double cavidad, double espesorAbs, int fraccion) {

        // Input variables del panel
        this.EspesorPanel = EspesorPanel / 1000;
        this.DistanciaRepeticion = DistanciaRepeticion / 1000;
        this.RadioPerf = RadioPerf / 1000;
        this.cavidad = cavidad / 1000;
        espesorAire = cavidad - espesorAbs;
        areaAbierta = (Math.PI * Math.pow(this.RadioPerf, 2)) / (Math.pow(this.DistanciaRepeticion, 2));
        endCorr = this.EspesorPanel + 2 * this.RadioPerf * (.8 * (1 - 1.47 * Math.sqrt(this.areaAbierta) + .47 * Math.sqrt(Math.pow(this.areaAbierta, 3))));

        // Vector frecuencia
        int largo = 0;
        double[] freq1 = {31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000};
        double[] freq3 = {31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000};

        switch (fraccion) {
            case 1:
                largo = 9;
                break;
            case 3:
                largo = 26;
                break;
        }

        this.freq = new double[largo];

        switch (fraccion) {
            case 1:
                System.arraycopy(freq1, 0, freq, 0, largo);
                break;
            case 3:
                System.arraycopy(freq3, 0, freq, 0, largo);
                break;
        }

        // Parametros y variables de ecuaciones
        ////  Doubles        
        this.w = new double[freq.length];
        double real[] = new double[freq.length];
        double img[] = new double[freq.length];
        double kaire[] = new double[freq.length];
        double capaviscosa[] = new double[freq.length];    //No se usa
        double capatermica[] = new double[freq.length];    //No se usa        
        double DyB[] = new double[freq.length];

        ////  Complejos
        Zc = new Complejo[freq.length];
        ZAire_c = new Complejo[freq.length];
        Kg = new Complejo[freq.length];
        CotKaireTaire = new Complejo[freq.length];
        CotKaTa = new Complejo[freq.length];
        J = new Complejo[freq.length];
        negJ = new Complejo[freq.length];
        Complejo CosKaTa[] = new Complejo[freq.length];
        Complejo SinKaTa[] = new Complejo[freq.length];

        // Calculos por banda de fcia 'i'
        for (int i = 0; i < freq.length; i++) {
            // Vector J
            J[i] = new Complejo(0, 1);
            negJ[i] = new Complejo(0, -1);

            w[i] = 2 * Math.PI * freq[i];

            kaire[i] = F.getTermInter1() * freq[i];

            //        Valores capa limite
            capaviscosa[i] = 1000 * Math.sqrt((2 * F.getViscoCinem()) / (F.getDensAire() * w[i]));
            capatermica[i] = 1000 * Math.sqrt((2 * F.getConduTerm()) / (F.getDensAire() * F.getCapCalorEspe() * w[i]));

            //        Valores absorbente
            DyB[i] = (F.getDensAire() * freq[i]) / F.getResisflujo();
            real[i] = F.getZAire() * (1 + 0.0571 * (Math.pow(DyB[i], -0.754)));
            img[i] = F.getZAire() * (-0.087 * (Math.pow(DyB[i], -0.732)));

            Zc[i] = new Complejo(real[i], img[i]);
            ZAire_c[i] = new Complejo();
            ZAire_c[i].setReal(F.getZAire());
            ZAire_c[i].setImag(0);

            real[i] = kaire[i] * (1 + 0.0978 * Math.pow(DyB[i], -.7));
            img[i] = kaire[i] * (-.189 * Math.pow(DyB[i], -.595));
            Kg[i] = new Complejo(real[i], img[i]);

            // Valores intermedios
            real[i] = Math.cos(kaire[i] * espesorAire / 1000) / (Math.sin(kaire[i] * espesorAire / 1000));
            CotKaireTaire[i] = new Complejo(real[i], 0);
            CosKaTa[i] = (Kg[i].multiplicar(espesorAbs / 1000)).coseno();
            SinKaTa[i] = (Kg[i].multiplicar(espesorAbs / 1000)).seno();
            CotKaTa[i] = CosKaTa[i].dividir(SinKaTa[i]);
        }
    }

    public double[] getFreq() {
        return freq;
    }
    
    public double getAreaAbierta() {
        return areaAbierta;
    }

    public double[] abs1() {

        double alpha[] = new double[freq.length];
        Complejo Z1[] = new Complejo[freq.length];
        Complejo Z2[] = new Complejo[freq.length];
        Complejo Z3[] = new Complejo[freq.length];
        Complejo R[] = new Complejo[freq.length];
        Complejo Zm[] = new Complejo[freq.length];

        // Cálculo impedancias
        // Absorbente contra Panel
        for (int i = 0; i < freq.length; i++) {

            Z1[i] = CotKaireTaire[i].multiplicar(F.getZAire()).multiplicar(negJ[i]);

            Z2[i] = ((negJ[i].multiplicar(Z1[i].multiplicar(Zc[i].multiplicar(CotKaTa[i]))))
                    .sumar(Zc[i].multiplicar(Zc[i]))).
                    dividir((Z1[i].restar((Zc[i].multiplicar(CotKaTa[i]).multiplicar(J[i])))));

            Zm[i] = new Complejo((F.getDensAire() / areaAbierta) * Math.sqrt(8 * F.getViscoCinem() * w[i])
                    * (1 + (endCorr / (2 * RadioPerf))), 0);

            Z3[i] = Z2[i].sumar(Zm[i]).sumar((J[i].dividir(areaAbierta)).multiplicar(endCorr * w[i] * F.getDensAire()));

            //  Cálculo coeficientes
            R[i] = (Z3[i].restar(ZAire_c[i])).dividir((Z3[i].sumar(ZAire_c[i])));
            alpha[i] = Math.round((1 - Math.pow(R[i].modulo(), 2)) * 100);
            alpha[i] /= 100;

            // Puede darse que algun alpha de negativo (ej, -0.01), por lo tanto convertimos estos posibles
            // valores en cero
            
            if (alpha[i] < 0) {
                alpha[i] = 0;
            }
        }
        return alpha;
    }

    public double[] abs2() {

        double alpha[] = new double[freq.length];
        Complejo Z1[] = new Complejo[freq.length];
        Complejo Z2[] = new Complejo[freq.length];
        Complejo Z3[] = new Complejo[freq.length];
        Complejo R[] = new Complejo[freq.length];

        // Cálculo impedancias
        // Absorbente contra Muro
        for (int i = 0; i < freq.length; i++) {

            Z1[i] = CotKaTa[i].multiplicar(Zc[i]).multiplicar(negJ[i]);

            Z2[i] = ((negJ[i].multiplicar(Z1[i].multiplicar(ZAire_c[i].multiplicar(CotKaireTaire[i]))))
                    .sumar(ZAire_c[i].multiplicar(ZAire_c[i]))).
                    dividir((Z1[i].restar((ZAire_c[i].multiplicar(CotKaireTaire[i]).multiplicar(J[i])))));

            Z3[i] = Z2[i].sumar(J[i].multiplicar(endCorr / areaAbierta * F.getDensAire() * w[i]));
            Z3[i].setReal(Z3[i].getReal() + F.getDensAire() / areaAbierta * Math.sqrt(8 * F.getViscoCinem() * w[i]) * ((EspesorPanel / 2 * RadioPerf) + 1));

            //  Cálculo coeficientes
            R[i] = (Z3[i].restar(ZAire_c[i])).dividir((Z3[i].sumar(ZAire_c[i])));
            alpha[i] = Math.round((1 - Math.pow(R[i].modulo(), 2)) * 100);
            alpha[i] /= 100;

            // Puede darse que algun alpha de negativo (ej, -0.01), por lo tanto convertimos estos posibles
            // valores en cero
            if (alpha[i] < 0) {
                alpha[i] = 0;
            }
        }
        return alpha;
    }

    public double[] abs3() {

        double alpha[] = new double[freq.length];
        Complejo Z1[] = new Complejo[freq.length];
        Complejo Z2[] = new Complejo[freq.length];
        Complejo CotKaD[] = new Complejo[freq.length];
        Complejo R[] = new Complejo[freq.length];

        // Cálculo impedancia
        // Sin cámara de aire
        for (int i = 0; i < freq.length; i++) {

            CotKaD[i] = ((Kg[i].multiplicar(cavidad)).coseno()).dividir(((Kg[i].multiplicar(cavidad)).seno()));

            Z1[i] = CotKaD[i].multiplicar(Zc[i]).multiplicar(negJ[i]);

            Z2[i] = ((J[i].multiplicar(w[i] * F.getDensAire() * endCorr / areaAbierta)
                    .sumar(Z1[i])));

            //  Cálculo coeficientes
            R[i] = (Z2[i].restar(ZAire_c[i])).dividir((Z2[i].sumar(ZAire_c[i])));
            alpha[i] = Math.round((1 - Math.pow(R[i].modulo(), 2)) * 100);
            alpha[i] /= 100;

            // Puede darse que algun alpha de negativo (ej, -0.01), por lo tanto convertimos estos posibles
            // valores en cero
            if (alpha[i] < 0) {
                alpha[i] = 0;
            }
        }
        return alpha;
    }
}
