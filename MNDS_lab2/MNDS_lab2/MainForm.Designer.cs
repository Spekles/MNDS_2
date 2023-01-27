
namespace MNDS_lab2
{
    partial class MainForm
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.outputTextBox = new System.Windows.Forms.RichTextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.firstDiameterTextBox = new System.Windows.Forms.TextBox();
            this.secondDiameterTextBox = new System.Windows.Forms.TextBox();
            this.thicknessTextBox = new System.Windows.Forms.TextBox();
            this.calculateButton = new System.Windows.Forms.Button();
            this.label4 = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // outputTextBox
            // 
            this.outputTextBox.Location = new System.Drawing.Point(12, 12);
            this.outputTextBox.Name = "outputTextBox";
            this.outputTextBox.Size = new System.Drawing.Size(420, 346);
            this.outputTextBox.TabIndex = 0;
            this.outputTextBox.Text = "";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(453, 31);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(193, 20);
            this.label1.TabIndex = 1;
            this.label1.Text = "Enter first diameter, mm, d0";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(454, 84);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(209, 20);
            this.label2.TabIndex = 2;
            this.label2.Text = "Enter second diameter, mm, D";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(453, 138);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(200, 20);
            this.label3.TabIndex = 3;
            this.label3.Text = "Enter thickness of ring, mm, s";
            // 
            // firstDiameterTextBox
            // 
            this.firstDiameterTextBox.Location = new System.Drawing.Point(454, 54);
            this.firstDiameterTextBox.Name = "firstDiameterTextBox";
            this.firstDiameterTextBox.Size = new System.Drawing.Size(125, 27);
            this.firstDiameterTextBox.TabIndex = 4;
            // 
            // secondDiameterTextBox
            // 
            this.secondDiameterTextBox.Location = new System.Drawing.Point(453, 108);
            this.secondDiameterTextBox.Name = "secondDiameterTextBox";
            this.secondDiameterTextBox.Size = new System.Drawing.Size(125, 27);
            this.secondDiameterTextBox.TabIndex = 5;
            // 
            // thicknessTextBox
            // 
            this.thicknessTextBox.Location = new System.Drawing.Point(454, 162);
            this.thicknessTextBox.Name = "thicknessTextBox";
            this.thicknessTextBox.Size = new System.Drawing.Size(125, 27);
            this.thicknessTextBox.TabIndex = 6;
            // 
            // calculateButton
            // 
            this.calculateButton.Location = new System.Drawing.Point(454, 196);
            this.calculateButton.Name = "calculateButton";
            this.calculateButton.Size = new System.Drawing.Size(94, 29);
            this.calculateButton.TabIndex = 7;
            this.calculateButton.Text = "Calculate";
            this.calculateButton.UseVisualStyleBackColor = true;
            this.calculateButton.Click += new System.EventHandler(this.calculateButton_Click);
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(454, 258);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(178, 60);
            this.label4.TabIndex = 8;
            this.label4.Text = "The work was done by\r\nVolodymyr Shyrokoradiuk\r\nKNm-22-1";
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 20F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(735, 376);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.calculateButton);
            this.Controls.Add(this.thicknessTextBox);
            this.Controls.Add(this.secondDiameterTextBox);
            this.Controls.Add(this.firstDiameterTextBox);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.outputTextBox);
            this.Name = "MainForm";
            this.Text = "Labaratory Work №2";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.RichTextBox outputTextBox;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox firstDiameterTextBox;
        private System.Windows.Forms.TextBox secondDiameterTextBox;
        private System.Windows.Forms.TextBox thicknessTextBox;
        private System.Windows.Forms.Button calculateButton;
        private System.Windows.Forms.Label label4;
    }
}

