using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MNDS_lab2
{
    class FEM
    {
        double pi = Math.PI;       
        double r0, R, rj, Rj;
        double nu = 0.3,

        E = 2.1E5,
        sig_ys = 160,
        sig_u = 280,
        e_u = 0.32,
        e_tot = 0.4;

        double du0 = 0.01;
        double[] u = new double[8];
        double[] sig_r = new double[8];
        double[] sig_t = new double[8];
        double[] sig_ = new double[8];
        double[] sig_f = new double[8];
        double[] eps_pl = new double[8];

        int j_step = 1, N_dr = 6;
        int N_step = 2;
        //Введення вихідних даних для методу скінченних елементів:
        int N_elem = 144; //Число елементів.
        int N_nodes = 91; //Число вузлів.

        //Аналітичний розв‘язок пружної задачі:
        public string Elasticity(double d0, double D, double s)
        {            
            string result;
            double dr, Cel, Cel_pl, du0_el, K0, Fmax;
            double[] r = new double[8];
            double[] eps_el_r = new double[8];
            double[] eps_el_t = new double[8];
            int i;
            du0_el = -du0;
            r0 = d0 / 2;
            R = D / 2;
            Rj = r0;
            Rj = R;
            dr = (R - r0) / N_dr;
            for (i = 1; i <= N_dr + 1; i++) //Організація циклу по і:
            {
                r[i] = r0 + (i - 1) * dr;
                K0 = D / d0;
                Cel = Math.Pow(K0, 2) * (1 + nu) / (1 - nu + Math.Pow(K0, 2) * (1 + nu));
                //Розрахунок переміщень:
                u[i] = du0_el / r0 * (r[i] - Cel * (r[i] - Math.Pow(r0, 2) / r[i]));
                //Деформацій:
                eps_el_r[i] = du0_el / r0 * (1 - Cel * (1 + Math.Pow(r0, 2) / Math.Pow(r[i], 2)));
                eps_el_t[i] = du0_el / r0 * (1 - Cel * (1 - Math.Pow(r0, 2) / Math.Pow(r[i], 2)));
                //Напруження:
                sig_r[i] = E / (1 - Math.Pow(nu, 2)) * (eps_el_r[i] + nu * eps_el_t[i]);
                sig_t[i] = E / (1 - Math.Pow(nu, 2)) * (eps_el_t[i] + nu * eps_el_r[i]);
                //Інтенсивність напружень:
                sig_[i] = Math.Sqrt(Math.Pow(sig_r[i], 2) - sig_r[i] * sig_t[i] + Math.Pow(sig_t[i], 2));
            }

            result = "u[1] (mm) = " + u[1].ToString() + "\n" + "u[7] (mm) = " + u[7].ToString() + "\n" + "\n";
            return result;

        }
        //Кінец функції Elasticity() для аналітичного розв‘язку.
        //Чисельний розв‘язок пружної задачі.
        public string Elasticity_FEM(double d0, double D, double s)
        {
            string result;
            double du0_el, rj, Rj;
            int i, j, l, m, n;
            du0_el = -du0;
            rj = r0; Rj = R;
            //- до 1-го етапу навантаження j_step=l.
            //Координати 13 вузлів на внутрішньому x_rj[t], y_rj[t](t-1,2,3,.-..,13) та
            //13 вузлів на зовнішньому x_Rj[t], y_Rj[t] контурах кільця

            double d_theta;
            double[] theta = new double[14];
            double[] x_rj = new double[14];
            double[] y_rj = new double[14];
            double[] x_Rj = new double[14];
            double[] y_Rj = new double[14];
            int t, N_d_theta;
            //Кількість секторів кута theta=pi/2 для 1 четвертини:
            N_d_theta = 12;
            //Для одного сектора кут дорівнює:
            d_theta = pi / 2 / N_d_theta;

            for (t = 1; t <= N_d_theta + 1; t++)
            {
                theta[t] = pi / 2 - (t - 1) * d_theta;
                x_rj[t] = rj * Math.Cos(theta[t]);
                y_rj[t] = rj * Math.Sin(theta[t]);
                x_Rj[t] = Rj * Math.Cos(theta[t]);
                y_Rj[t] = Rj * Math.Sin(theta[t]);
            }
            //Кінець розрахунку координат вузлів на контурах кільця
            //Довжина d_jiorm(t] шести відрізків t-x нормалей,
            //що виходять з t-x вузлів внутрішнього контуру:

            double[] d_norm = new double[14];
            for (t = 1; t <= N_d_theta + 1; t++)
                d_norm[t] = Math.Sqrt((x_Rj[t] - x_rj[t]) * (x_Rj[t] - x_rj[t]) + (y_Rj[t] - y_rj[t]) * (y_Rj[t] - y_rj[t])) / 6;

            //Координати x_glob[g],y_glob[g] g-x вузлів
            //з глобальними номерами g=-1,2, 3, ..., N_nodes:
            double[] x_glob = new double[92];
            double[] y_glob = new double[92];
            int g, k1;
            i = -6; j = 0;
            for (t = 1; t <= N_d_theta + 1; t++)
            {
                k1 = 0; i = i + 7; j = j + 7;
                for (g = i; g <= j; g++)
                {
                    k1 = k1 + 1;
                    x_glob[g] = x_rj[t] + (k1 - 1) * d_norm[t] * Math.Cos(theta[t]);
                    y_glob[g] = y_rj[t] + (k1 - 1) * d_norm[t] * Math.Sin(theta[t]);
                }
            }
            //Завершення розрахунку координат вузлів з глобальними номерами.

            //Розрахунок матриці індексів.
            //Введення першої частини матриці індексів для першої
            //мінімальної та симетричної частини сітки скінченних елементів
            // (п=1,2,3, . .,24):
            int[,] mi = new int[145, 4];
            int[,] MI = new int[145, 7];
            mi[1, 1] = 1; mi[1, 2] = 8; mi[1, 3] = 9;
            mi[2, 1] = 1; mi[2, 2] = 9; mi[2, 3] = 2;
            mi[3, 1] = 2; mi[3, 2] = 9; mi[3, 3] = 3;
            mi[4, 1] = 3; mi[4, 2] = 9; mi[4, 3] = 10;
            mi[5, 1] = 3; mi[5, 2] = 10; mi[5, 3] = 11;
            mi[6, 1] = 3; mi[6, 2] = 11; mi[6, 3] = 4;
            mi[7, 1] = 4; mi[7, 2] = 11; mi[7, 3] = 5;
            mi[8, 1] = 5; mi[8, 2] = 11; mi[8, 3] = 12;
            mi[9, 1] = 5; mi[9, 2] = 12; mi[9, 3] = 13;
            mi[10, 1] = 5; mi[10, 2] = 13; mi[10, 3] = 6;
            mi[11, 1] = 6; mi[11, 2] = 13; mi[11, 3] = 7;
            mi[12, 1] = 7; mi[12, 2] = 13; mi[12, 3] = 14;
            mi[13, 1] = 8; mi[13, 2] = 15; mi[13, 3] = 9;
            mi[14, 1] = 9; mi[14, 2] = 15; mi[14, 3] = 16;
            mi[15, 1] = 9; mi[15, 2] = 16; mi[15, 3] = 17;
            mi[16, 1] = 9; mi[16, 2] = 17; mi[16, 3] = 10;
            mi[17, 1] = 10; mi[17, 2] = 17; mi[17, 3] = 11;
            mi[18, 1] = 11; mi[18, 2] = 17; mi[18, 3] = 18;
            mi[19, 1] = 11; mi[19, 2] = 18; mi[19, 3] = 19;
            mi[20, 1] = 11; mi[20, 2] = 19; mi[20, 3] = 12;
            mi[21, 1] = 12; mi[21, 2] = 19; mi[21, 3] = 13;
            mi[22, 1] = 13; mi[22, 2] = 19; mi[22, 3] = 20;
            mi[23, 1] = 13; mi[23, 2] = 20; mi[23, 3] = 21;
            mi[24, 1] = 13; mi[24, 2] = 21; mi[24, 3] = 14;

            //Розрахунок іншої частини матриці індексів
            //для іншої частини сітки скінченних елементів (п=25.. .144),
            //використовуючи першу частину матриці індексів для 24 скінченних
            //елементів:
            int n1, n2;
            for (i = 1; i <= 5; i++)
            {
                n1 = 1 + i * 24; n2 = n1 + 23;
                for (n = n1; n <= n2; n++)
                {
                    for (l = 1; l <= 3; l++)
                    {
                        mi[n, l] = mi[n - 24, l] + 14;
                    }
                }
            }
            //Завершення розрахунку матриці індексів mi

            //Матриця індексів МІ для степенів вільності:
            for (n = 1; n <= N_elem; n++)
            {
                l = 1;
                MI[n, 1] = mi[n, l] * 2 - 1;
                MI[n, 2] = mi[n, l] * 2;

                l = 2;
                MI[n, 3] = mi[n, l] * 2 - 1;
                MI[n, 4] = mi[n, l] * 2;

                l = 3;
                MI[n, 5] = mi[n, l] * 2 - 1;
                MI[n, 6] = mi[n, l] * 2;
            }

            //координати вузлів x[n][l], y[n][l]
            //з локальними номерами l=1,2,3 для n–х елементів:
            //n=1,2,3,…,N_elem:
            double[,] x = new double[145, 4];
            double[,] y = new double[145, 4];
            for (n = 1; n <= N_elem; n++)
            {
                for (l = 1; l <= 3; l++)
                {
                    g = mi[n, l];
                    x[n, l] = x_glob[g];
                    y[n, l] = y_glob[g];
                }
            }
            //Завершення розрахунку координат вузлів з локальними номерами

            //Матриця пружності
            double[,] Dp = new double[4, 4];
            for (i = 1; i <= 3; i++)
                for (j = 1; j <= 3; j++)
                    Dp[i, j] = 0;
            Dp[1, 1] = Dp[2, 2] = 1;
            Dp[2, 1] = Dp[1, 2] = nu;
            Dp[3, 3] = (1 - nu) / 2;
            for (i = 1; i <= 3; i++)
                for (j = 1; j <= 3; j++)
                    Dp[i, j] = Dp[i, j] * E / (1 - Math.Pow(nu, 2));
            //Завершення розрахунку матриці пружності

            //Розрахунок локальних k та глобальних матриць жорсткості
            double A_elem;
            double[,] k = new double[7, 7];
            double[,] K = new double[183, 184];
            double[,] B = new double[4, 7];
            double[,] B_T = new double[7, 4];
            double[,] B_T_D = new double[7, 4];
            double[] s_elem = new double[145];
            int i_glob, j_glob;

            //Обнуління глобальної матриці жорсткості
            //(стовпець "+1" - для вектора вузлових сил):
            for (i = 1; i <= 2 * N_nodes; i++)
                for (j = 1; j <= N_nodes; j++)
                    K[i, j] = 0;

            //обнуління матриці диференційних переміщень:
            for (i = 1; i <= 3; i++)
                for (j = 1; j <= 6; j++)
                    B[i, j] = 0;

            //Задача товщини скінченних елементів:
            for (n = 1; n <= N_elem; n++)
                s_elem[n] = s;

            //Основний цикл no n для розрахунку глобальної матриці жорсткості:
            for (n = 1; n <= N_elem; n++)
            {
                //Обнуління локальної матриці жорсткості:
                for (i = 1; i <= 6; i++)
                    for (j = 1; j <= 6; j++)
                        k[i, j] = 0;

                //Розрахунок матриці диференціальних переміщень:
                B[1, 1] = B[3, 2] = y[n, 2] - y[n, 3];
                B[1, 3] = B[3, 4] = y[n, 3] - y[n, 1];
                B[1, 5] = B[3, 6] = y[n, 1] - y[n, 2];
                B[2, 2] = B[3, 1] = x[n, 3] - x[n, 2];
                B[2, 4] = B[3, 3] = x[n, 1] - x[n, 3];
                B[2, 6] = B[3, 5] = x[n, 2] - x[n, 1];

                //Площа n-гo скінченного елемента:
                A_elem = (x[n, 2] * y[n, 3] + x[n, 1] * y[n, 2] + y[n, 1] * x[n, 3] - y[n, 1] * x[n, 2] - y[n, 2] * x[n, 3] - x[n, 1] * y[n, 3]) / 2;
                for (i = 1; i <= 3; i++)
                    for (j = 1; j <= 6; j++)
                        B[i, j] = B[i, j] / (2 * A_elem);
                //Кінець В

                //Транспонована матриця В_Т:
                for (i = 1; i <= 6; i++)
                    for (j = 1; j <= 3; j++)
                        B_T[i, j] = B[j, i];

                //Множення матриць В__Т и D, та отримання В_Т_D:
                for (i = 1; i <= 6; i++)
                    for (j = 1; j <= 3; j++)
                        B_T_D[i, j] = 0;

                for (i = 1; i <= 6; i++)
                    for (j = 1; j <= 3; j++)
                        for (m = 1; m <= 3; m++)
                            B_T_D[i, j] = B_T_D[i, j] + B_T[i, m] * Dp[m, j];

                //Локальна матриця жорсткості k для n-го елемента:
                for (i = 1; i <= 6; i++)
                    for (j = 1; j <= 6; j++)
                        for (m = 1; m <= 3; m++)
                            k[i, j] = k[i, j] + B_T_D[i, m] * B[m, j];

                for (i = 1; i <= 6; i++)
                    for (j = 1; j <= 6; j++)
                        k[i, j] = k[i, j] * A_elem * s_elem[n];
                //Кінець k.

                //Сборка глобальної матриці жорсткості:
                for (i = 1; i <= 6; i++)
                {
                    i_glob = MI[n, i];
                    for (j = 0; j <= 6; j++)
                    {
                        j_glob = MI[n, j];                        
                        K[i_glob, j_glob] = K[i_glob, j_glob] + k[i, j];
                    }
                    //The end for j.
                }
                //The end for i.
            }
            //Звершення основного циклу по п для розрахунку
            //глобальної матриці жорсткості

            //Зборка глобального вектора вузлових сил F[183]
            //та урахування граничних умов:
            double[] F = new double[183];
            for (i = 1; i <= 2 * N_nodes; i++)
                F[i] = 0;

            //Урахування заданих ненульових переміщень du0_el
            //на внутрішньому контурі кільця:
            du0_el = -du0; t = 0;
            for (j = 1; j <= 169; j += 14)
            {
                t = t + 1;
                for (i = 1; i <= 2 * N_nodes; i++)
                {
                    F[i] = F[i] - K[i, j] * du0_el * Math.Cos(theta[t]);
                }
            }

            t = 0;
            for (j = 2; j <= 170; j += 14)
            {
                t = t + 1;
                for (i = 1; i <= 2 * N_nodes; i++)
                {
                    F[i] = F[i] - K[i, j] * du0_el * Math.Sin(theta[t]);
                }
            }

            //Визначення переміщень у векторі сил:
            t = 0;
            for (j = 1; j <= 169; j += 14)
            {
                t = t + 1;
                F[j] = du0_el * Math.Cos(theta[t]);
            }
            t = 0;
            for (j = 2; j <= 170; j += 14)
            {
                t = t + 1;
                F[j] = du0_el * Math.Sin(theta[t]);
            }
            //Завершення урахування заданих ненульових переміщень du0_l.

            //Урахування нульових тангенціальних переміщень на осі у:
            for (j = 1; j <= 13; j += 2) F[j] = 0;

            //Урахування нульових тангенціальних переміщень на осі х:
            for (j = 170; j <= 182; j += 2) F[j] = 0;
            //Завершення розрахунку глобального вектора вузлових сил.

            //Перетворення глобальної матриці жорсткості
            //з урахуванням граничних умов
            //На внутрішньому контурі кільця:
            for (j = 1; j <= 169; j += 14)
            {
                for (i = 1; i <= 2 * N_nodes; i++)
                {
                    if (i == j)
                    {
                        K[i, j] = 1;
                    }
                    else
                    {
                        K[i, j] = K[j, i] = 0;
                    }
                }
            }

            for (j = 2; j <= 17; j += 14)
            {
                for (i = 1; i <= 2 * N_nodes; i++)
                {
                    if (i == j)
                    {
                        K[i, j] = 1;
                    }
                    else
                    {
                        K[i, j] = 0;
                    }
                }
            }

            for (j = 3; j <= 13; j += 2) //На осі у:
            {
                for (i = 1; i <= 2 * N_nodes; i++)
                {
                    if (i == j)
                    {
                        K[i, j] = 1;
                    }
                    else
                    {
                        K[i, j] = K[j, i] = 0;
                    }
                }
            }

            for (j = 172; j <= 182; j += 2) //На осiи x:
            {
                for (i = 1; i <= 2 * N_nodes; i++)
                {
                    if (i == j)
                    {
                        K[i, j] = 1;
                    }
                    else
                    {
                        K[i, j] = K[j, i] = 0;
                    }
                }
            }
            //Завершення перетворення глобальної матриці жорсткості
            //з урахуванням граничних умов.

            //Визначення вектора сил F у вигляді останнього стовпця
            //глобальної матриці жорсткості
            j = 2 * N_nodes + 1;
            for (i = 1; i <= 2 * N_nodes; i++) K[i, j] = F[i];

            //Розв’язок системи рівнянь методом Гауса
            //та розрахунок глобального вектора переміщень U1183):
            double temp;
            int u1, N;
            double[] U = new double[183];
            u1 = 0; N = 2 * N_nodes;
            
        M0:
            u1 = u1 + 1;

            for (k1 = u1; k1 <= N; k1++)
            {
                if (K[k1, u1] != 0) goto M1;
            }

            result = "Error - zero";
            return result;
        M1:
            if (k1 == u1) goto M2;

            for (m = u1; m <= N + 1; m++)
            {
                temp = K[u1, m];
                K[u1, m] = K[k1, m];
                K[k1, m] = temp;
            }
        M2:
            for (j = N + 1; j >= u1; j--)
                K[u1, j] = K[u1, j] / K[u1, u1];

            for (i = k1 + 1; i <= N; i++)
            {
                for (j = u1 + 1; j <= N + 1; j++)
                {
                    K[i, j] = K[i, j] - K[i, u1] * K[u1, j];
                }
            }

            if (u1 != N) goto M0;

            for (i = N; i >= 1; i--)
            {
                U[i] = K[i, N + 1] / K[i, i];
                for (k1 = i - 1; k1 >= 1; k1--)
                {
                    K[k1, N + 1] = K[k1, N + 1] - K[k1, i] * U[i];
                }
            }
        Exit: //Завершення розв‘язку системи рівнянь методом Гауса
              //Перехід від глобальних U[183] до локальних u[(п][l]
              //переміщень для n-го скінченного елемента:
            double[,] u = new double[145, 7];
            for (n = 1; n <= N_elem; n++)
            {
                for (l = 1; l <= 6; l++)
                {
                    i_glob = MI[n, l];
                    U[l] = U[i_glob];
                }
            }

            //Об‘явлення матриці деформацій:
            double[,] eps_el = new double[145, 4];

            //Обнуління матриці деформацій:
            for (n = 1; n <= N_elem; n++)
                for (i = 1; i <= 3; i++)
                    eps_el[n, i] = 0;

            //Обнуління матриці диференціонування переміщень:
            for (i = 1; i <= 3; i++)
                for (j = 1; j <= 6; j++) B[i, j] = 0;

            //Розрахунок деформацій:
            for (n = 1; n <= N_elem; n++)
            {
                //Розрахунок матриці диференціонування переміщень:
                B[1, 1] = B[3, 2] = y[n, 2] - y[n, 3];
                B[1, 3] = B[3, 4] = y[n, 3] - y[n, 1];
                B[1, 5] = B[3, 6] = y[n, 1] - y[n, 2];
                B[2, 2] = B[3, 1] = x[n, 3] - x[n, 2];
                B[2, 4] = B[3, 3] = x[n, 1] - x[n, 3];
                B[2, 6] = B[3, 5] = x[n, 2] - x[n, 1];

                //Площа n-гo скінченного елемента:
                A_elem = (x[n, 2] * y[n, 3] + x[n, 1] * y[n, 2] + y[n, 1] * x[n, 3] - y[n, 1] * x[n, 2] - y[n, 2] * x[n, 3] - x[n, 1] * y[n, 3]) / 2;

                for (i = 1; i <= 3; i++)
                    for (j = 1; j <= 6; j++)
                        B[i, j] = B[i, j] / (2 * A_elem);
                //Завершення В.

                for (j = 1; j <= 3; j++)
                    for (i = 1; i <= 6; i++)
                        eps_el[n, j] = eps_el[n, j] + B[j, i] * u[n, i];
            }

            //Запис деформацій в координатах х,у:
            double[] eps_x = new double[145];
            double[] eps_y = new double[145];
            double[] gamma_xy = new double[145];

            for (n = 1; n <= N_elem; n++)
            {
                eps_x[n] = eps_el[n, 1];
                eps_y[n] = eps_el[n, 2];
                gamma_xy[n] = eps_el[n, 3];
            }

            //Об‘явлення матриці напружень:
            double[,] sig = new double[145, 4];

            //Обнуління матриці напружень:
            for (n = 1; n <= N_elem; n++)
                for (i = 1; i <= 3; i++)
                    sig[n, i] = 0;

            //Розрахунок напружень:
            int ll;
            for (n = 1; n <= N_elem; n++)
                for (k1 = 1; k1 <= 3; k1++)
                    for (ll = 1; ll <= 3; ll++)
                        sig[n, k1] = sig[n, k1] + Dp[k1, ll] * eps_el[n, ll];

            //Запис напружень в координатах х, у:
            double[] sig_x = new double[145];
            double[] sig_y = new double[145];
            double[] tau_xy = new double[145];

            for (n = 1; n <= N_elem; n++)
            {
                sig_x[n] = sig[n, 1];
                sig_y[n] = sig[n, 2];
                tau_xy[n] = sig[n, 3];
            }
            //Завершення розрахунку напружень в координатах х,у.

            //Розрахунок середніх деформацій та напружень у вузлах
            // в прямокутних координатах х,у.
            //Число скінченних елементів N_elem_node[g],
            //що сходяться в кожному вузлі
            // з глобальним номером g=1,2,3, ...,N_nodes:
            double[] N_elem_node = new double[92];
            int e;
            for (g = 1; g <= N_nodes; g++)
            {
                e = 0;
                for (n = 1; n <= N_elem; n++)
                {
                    for (l = 1; l <= 3; l++)
                    {
                        if (mi[n, 1] == g)
                        {
                            E = e + l; N_elem_node[g] = e;
                        }
                    }
                }
            }

            //Сумарні (sum) деформації та напруження у вузлі
            //з глобальним номером g=1,2.3,...,N_nodes:
            double[] eps_x_sum = new double[92];
            double[] eps_y_sum = new double[92];
            double[] gamma_xy_sum = new double[92];
            double[] sig_x_sum = new double[92];
            double[] sig_y_sum = new double[92];
            double[] tau_xy_sum = new double[92];

            //Обнуління масивів:
            for (g = 1; g <= N_nodes; g++)
                eps_x_sum[g] = eps_y_sum[g] = gamma_xy_sum[g] = sig_x_sum[g] = sig_y_sum[g] = tau_xy_sum[g] = 0;

            //Сумарні (sum) деформації та напруження:
            for (g = 1; g <= N_nodes; g++)
            {
                for (n = 1; n <= N_elem; n++)
                {
                    for (l = 1; l <= 3; l++)
                    {

                        if (mi[n, l] == g)
                        {
                            eps_x_sum[g] = eps_x_sum[g] + eps_x[n];
                            eps_y_sum[g] = eps_y_sum[g] + eps_y[n];
                            gamma_xy_sum[g] = gamma_xy_sum[g] + gamma_xy[n];
                            sig_x_sum[g] = sig_x_sum[g] + sig_x[n];
                            sig_y_sum[g] = sig_y_sum[g] + sig_y[n];
                            tau_xy_sum[g] = tau_xy_sum[g] + tau_xy[n];
                        }
                    }
                }
            }

            //Середі (mean) деформації та напруження у вузлах в прямокутних координатах ху:
            double[] eps_x_mean = new double[92];
            double[] eps_y_mean = new double[92];
            double[] gamma_xy_mean = new double[92];
            double[] sig_x_mean = new double[92];
            double[] sig_y_mean = new double[92];
            double[] tau_xy_mean = new double[92];

            for (g = 1; g <= N_nodes; g++)
            {
                eps_x_mean[g] = eps_x_sum[g] / N_elem_node[g];
                eps_y_mean[g] = eps_y_sum[g] / N_elem_node[g];
                gamma_xy_mean[g] = gamma_xy_sum[g] / N_elem_node[g];
                sig_x_mean[g] = sig_x_sum[g] / N_elem_node[g];
                sig_y_mean[g] = sig_y_sum[g] / N_elem_node[g];
                tau_xy_mean[g] = tau_xy_sum[g] / N_elem_node[g];
            }

            //Завершення розрахунку середніх деформацій та напружень
            //в координатах х, у.
            //Розрахунок середніх деформацій та напружень у вузлах
            //в полярних координатах r,theta.
            //Розрахунок полярного кута theta_glob[g] в g-му вузлі:
            double[] theta_glob = new double[192];
            i = -6; j = 0;
            for (t = 1; t <= N_d_theta + 1; t++)
            {
                i = i + 7; j = j + 7;
                for (g = i; g <= j; g++)
                {
                    theta_glob[g] = pi / 2 - (t - 1) * d_theta;
                }
            }

            //Середні (mean) деформації та напруження у вузлах
            //в полярних координатах r,theta:
            double[] eps_r_mean = new double[92];
            double[] eps_t_mean = new double[92];
            double[] gаmma_rt_mean = new double[92];
            double[] sig_r_mean = new double[92];
            double[] sig_t_mean = new double[92];
            double[] tau_rt_mean = new double[92];

            for (g = 1; g <= N_nodes; g++)
            {
                eps_r_mean[g] = (eps_x_mean[g] + eps_y_mean[g]) / 2 + (eps_x_mean[g] - eps_y_mean[g]) *
                                Math.Cos(2 * theta_glob[g]) / 2 + gamma_xy_mean[g] * Math.Sin(2 * theta_glob[g]);

                eps_t_mean[g] = (eps_x_mean[g] + eps_y_mean[g]) / 2 - (eps_x_mean[g] - eps_y_mean[g]) *
                                Math.Cos(2 * theta_glob[g]) / 2 - gamma_xy_mean[g] * Math.Sin(2 * theta_glob[g]);

                gаmma_rt_mean[g] = (eps_x_mean[g] - eps_y_mean[g]) * Math.Sin(2 * theta_glob[g]) / 2 *
                                    gamma_xy_mean[g] * Math.Cos(2 * theta_glob[g]);

                sig_r_mean[g] = (sig_x_mean[g] + sig_y_mean[g]) / 2 + (sig_x_mean[g] - sig_y_mean[g]) *
                                Math.Cos(2 * theta_glob[g]) / 2 + tau_xy_mean[g];

                sig_t_mean[g] = (sig_x_mean[g] + sig_y_mean[g]) / 2 - (sig_x_mean[g] - sig_y_mean[g]) *
                                Math.Cos(2 * theta_glob[g]) / 2 - tau_xy_mean[g] * Math.Sin(2 * theta_glob[g]);

                tau_rt_mean[g] = (sig_x_mean[g] - sig_y_mean[g]) * Math.Sin(2 * theta_glob[g]) / 2 + tau_xy_mean[g] *
                                 Math.Cos(2 * theta_glob[g]);
            }
            //Завершення розрахунку середніх деформацій та напружень
            //у вузлі в полярних координатах r,theta.

            //Інтенсивність деформацій:
            double[] eps__ = new double[92];
            for (g = 1; g <= N_nodes; g++)
                eps__[g] = 2 / Math.Sqrt(3) * Math.Sqrt(Math.Pow(eps_r_mean[g], 2) +
                eps_r_mean[g] * eps_t_mean[g] + Math.Pow(eps_t_mean[g], 2));

            //Інтенсивність напружень:
            double[] sig__ = new double[192];
            for (g = 1; g <= N_nodes; g++)
                sig__[g] = Math.Sqrt(Math.Pow(sig_r_mean[g], 2) -
                sig_r_mean[g] * sig_t_mean[g] + Math.Pow(sig_t_mean[g], 2));

            //Сила F;
            double dF, F1, F_def;
            double[] dl_rj = new double[13];
            i = -6; F1 = 0;
            for (t = 1; t <= N_d_theta; t++)
            {
                i = i + 7;
                dl_rj[t] = Math.Sqrt((x_rj[t + 1] - x_rj[t]) * (x_rj[t + 1] - x_rj[t]) +
                (y_rj[t + 1] - y_rj[t]) * (y_rj[t + 1] - y_rj[t]));
                dF = dl_rj[t] * s * (sig_r_mean[i] + sig_r_mean[i + 7]) / 2;
                F1 = F1 + dF;
            }
            F_def = 4 * F1;
            //Завершення розрахунку сили.

            //Переміщення u_Rj[t] на зовнішньому контурі кільця Rj:
            i = 14;
            double[] u_Rj = new double[14];
            for (t = 2; t <= 12; t++)
            {
                i = i + 14;
                u_Rj[t] = U[i - 1] * Math.Cos(theta[t]) + U[i] * Math.Sin(theta[t]) / 2;
            }
            u_Rj[1] = U[14]; u_Rj[13] = U[181];

            //Виведення результатів розрахунків
            //Виведення переміщень u_Rj[t]:
            result = "";
            for (t = 1; t <= N_d_theta + 1; t++)
            {
                result += "t = " + t.ToString() + ", u_rj[t] (mm) = " + u_Rj[t].ToString() + "\n";
            }

            return result; 
        }
        //завершення функції Elastlcity_FEM().      
    }
}
