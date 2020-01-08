/*
 * This file is part of THESIAS
 * Copyright (C) 2004-2020 David-Alexandre Trégouët, Valérie Garelle
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import javax.swing.JOptionPane;
import javax.swing.JDialog;
import javax.swing.JTextField;
import java.beans.*; 
import java.awt.*;
import java.awt.event.*;
import javax.swing.JLabel;
import java.io.*;
import javax.swing.*;
import javax.swing.filechooser.*;

class CustomDialog extends JFrame
                   implements ActionListener,
                              PropertyChangeListener {
    private String typedText = null;
    private JTextField textField;

	JFileChooser fc;
	private String magicWord;
	private JOptionPane optionPane;

	private String btnString1 = "Enter";
	private String btnString2 = "Cancel";
	private String sbtnClose = "Close";
	private String sbtnAdvance = "Advance";
	private String sbtnPrint   = "Print";
	private String sbtnSave    = "Save";

	JButton close;
	JButton advanceButton,printButton,saveButton;
	String sParam;
	Object covar[];
	GraficT grafic;
	dataThesias datta;
	private String sPrint;
	String sSave;
  


    public String getValidatedText() {
        return typedText;
    }

    /** Creates the reusable dialog. */
    public CustomDialog(dataThesias ddd,GraficT art,Frame aFrame, String aWord,StringBuffer s,String sParam,Object []covar,StringBuffer sbufTest) {
		super();
		sSave = new String(s);
		sPrint = new String(sbufTest);
		grafic = art;
		datta = ddd;
		this.covar = covar;

		fc = new JFileChooser();
		Color saumon = new Color( 255, 204, 153 );	
		setDefaultLookAndFeelDecorated(true);
		GridLayout grid = new GridLayout(2,0);

		getContentPane().setLayout(new BorderLayout());

		JEditorPane a = new JEditorPane();
		a.setContentType("text/html");
		a.setEditable(false);
		a.setText(s.toString());

		JScrollPane jscrol = new JScrollPane(a);
		getContentPane().add(jscrol,BorderLayout.CENTER);

		JPanel buttPanel = new JPanel(new GridLayout(1,2));
		close = new JButton(sbtnClose);
		close.setBackground(saumon);
		close.setActionCommand(sbtnClose);
		close.addActionListener(this);
					buttPanel.add(close);
		this.sParam = sParam;
		if (sParam.length() != 0)
		{
			
			advanceButton = new JButton("Advance");
			advanceButton.setBackground(saumon);
			advanceButton.setActionCommand(sbtnAdvance);
			advanceButton.addActionListener(this);
					buttPanel.add(advanceButton);

		}

		printButton = new JButton(sbtnPrint);
		printButton.setActionCommand(sbtnPrint);
		printButton.addActionListener(this);
		buttPanel.add(printButton);
		
		saveButton = new JButton(sbtnSave);
		saveButton.setActionCommand(sbtnSave);
		saveButton.addActionListener(this);
		buttPanel.add(saveButton);
		
		getContentPane().add(buttPanel,BorderLayout.SOUTH);
	
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
	    GraphicsDevice[] gs = ge.getScreenDevices();
    
	    // Get size of each screen
	
		DisplayMode dm = gs[0].getDisplayMode();
	    int screenWidth = dm.getWidth();
		int screenHeight = dm.getHeight();
	    
		this.setSize(dm.getWidth(),dm.getHeight()-80);

		printButton.setBackground(saumon);
		saveButton.setBackground(saumon);



    }

    public void actionPerformed(ActionEvent e) {

		if( sbtnClose.equals(e.getActionCommand()) ){
			setVisible(false);
		}
		else if (sbtnAdvance.equals(e.getActionCommand()) )
		{
			OptionDialog optDialog = new OptionDialog(datta,grafic,this,sParam,this.covar);
			optDialog.setVisible(true);
		}
		else if ( sbtnPrint.equals(e.getActionCommand()) )
		{
			TextPrint t = new TextPrint();
			t.runPrint(sPrint);
		}
		else if( sbtnSave.equals(e.getActionCommand()) )
		{
			int returnVal = fc.showSaveDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
				String content;
				try
				{
					int i,j,tj = 0; 
					FileWriter fwData;
					fwData = new FileWriter(file);
					fwData.write(sSave, 0 , sSave.length());
					fwData.close();
					fwData.close();
				}
				catch ( IOException eIO	)
				{
					System.out.println(eIO.toString());
				}
			}
		}
    }

    public void propertyChange(PropertyChangeEvent e) {
  
    }

   
    public void clearAndHide() {
        textField.setText(null);
        setVisible(false);
    }
}
