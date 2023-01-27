using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace MNDS_lab2
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
        }

        private void calculateButton_Click(object sender, EventArgs e)
        {
            string res1, res2;
            double d0, D, s;

            if (!Double.TryParse(firstDiameterTextBox.Text,out d0) ||
                !Double.TryParse(secondDiameterTextBox.Text, out D) ||
                !Double.TryParse(thicknessTextBox.Text, out s))
            {
                MessageBox.Show("Put correct data in program!", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            FEM fem = new FEM();
            outputTextBox.Text = "";

            res1 = fem.Elasticity(d0, D, s);
            res2 = fem.Elasticity_FEM(d0, D, s);
            outputTextBox.Text = res1 + res2;
        }
    }
}
