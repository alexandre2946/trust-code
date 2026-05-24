Example:

    Statistiques Dt_post dtst {

      t_deb 0.1 t_fin 0.12

      Moyenne Pression

      Ecart_type Pression

      Correlation Vitesse Vitesse

    }

will write every dt_post the mean, standard deviation and correlation value:

if \f$t < t_{deb}\f$ or \f$t > t_{fin}\f$

\f[ \text{average: } \overline{P(t)} = 0 \f]
\f[ \text{std deviation: } \langle P(t) \rangle = 0 \f]
\f[ \text{correlation: } \langle U(t) \cdot V(t) \rangle = 0 \f]

if \f$t > t_{deb}\f$ and \f$t < t_{fin}\f$

\f[ \text{average: } \overline{P(t)} = \frac{1}{t - t_{deb}} \int_{t_{deb}}^{t} P(s)\, ds \f]
\f[ \text{std deviation: } \langle P(t) \rangle = \sqrt{ \frac{1}{t - t_{deb}} \int_{t_{deb}}^{t} \left[ P(s) - \overline{P(t)} \right]^2 ds } \f]
\f[ \text{correlation: } \langle U(t) \cdot V(t) \rangle = \frac{1}{t - t_{deb}} \int_{t_{deb}}^{t} \left[ U(s) - \overline{U(t)} \right] \cdot \left[ V(s) - \overline{V(t)} \right] ds \f]
