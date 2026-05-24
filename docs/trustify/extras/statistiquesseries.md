Example:

    Statistiques_en_serie Dt_integr dtst {

      Moyenne Pression

    }

will calculate and write every dtst seconds the mean value:

\f[ (n+1)\,dt_{integr} > t > n \cdot dt_{integr}, \quad \overline{P(t)} = \frac{1}{t - n \cdot dt_{integr}} \int_{n \cdot dt_{integr}}^{t} P(t)\, dt \f]
