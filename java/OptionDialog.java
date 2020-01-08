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
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.ImageIcon;
import javax.swing.BoxLayout;
import javax.swing.Box;
import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.JTabbedPane;
import javax.swing.JPanel;
import javax.swing.JFrame;
import java.beans.*; 
import java.awt.*;
import java.awt.event.*;
import javax.swing.JTable;
import javax.swing.table.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.datatransfer.*;
import java.util.*;
import java.io.*;
import javax.swing.JComboBox;


class OptionDialog extends JFrame implements ActionListener, PropertyChangeListener {
	private String magicWord;
	private JOptionPane optionPane;	
	private String sbtnApply = "Run";
	private String sbtnCancel = "Cancel";
	private JButton jbtnApply,jbtnCancel;
	private String sHypothesis	   = "Hypothesis";
	private String sHypothesisTip  = "Do you want to test specific hypothesis on haplotypic effects";
	private String sHomogeneity    = "Homogeneity";
	private String sHomogeneityTip = "Do you want to test for the homogeneity of some allelic effects";
	private String sNonAdditivity = "Non-additivity";
	private String sNonAdditivityTip ="Do you want to test the non-additivity of some haplotypic effects";
	private String sEnvironmentInteractions = "Covariate interactions"; 
	private String sEnvironmentInteractionsTip = "Do you want to test haplotype * environment interactions ";
	private Frame parent;
	Hipothesis	hypothesis;
	private Hipothesis aditivity;
	Hipothesis	interaction;
	Homogeneity homogeneity;
	GraficT grafic;
	dataThesias datta;
	private Cursor hourglassCursor;
	private Cursor normalCursor;
    public OptionDialog(dataThesias datta,GraficT art,Frame aFrame,String sParam,Object []covar) {
		super();
		Color chair = new Color( 249, 230, 196 );
		Color saumon = new Color( 255, 204, 153 );
		super.setBackground(chair);
		hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
		normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
		parent = aFrame;
		this.datta = datta;
		grafic = art;
		int nkat = 1;
		if ( grafic.kindOfPhenotype == 5 )
		{
			String content_prec = null;
			String content = null;
			try
			{
				FileReader s = new FileReader("para.txt");
				BufferedReader b = new BufferedReader(s);
	
				while( (content = b.readLine())!= null ){content_prec = content;}
				b.close();
				s.close();
			}
			catch ( IOException eIO	){System.out.println(eIO.toString());}
			nkat = Integer.parseInt(content_prec);
			System.out.println("----------------------------");
			System.out.println("nkat = "+nkat);
			System.out.println("sParam : " +sParam);
		}
		setBounds(10,10,100,30);
		StringTokenizer sTok = new StringTokenizer(sParam);
		Object dataTemp[] = new Object[1000];
		Object dataBun[];
		if ( sTok.hasMoreTokens() )
			sTok.nextToken();
		int i = 0, k =0;
		String a,b;
		while( sTok.hasMoreTokens() )
		{
			i++;
			sTok.nextToken();
			sTok.nextToken();
			a = sTok.nextToken();
			b = sTok.nextToken();
			
			if( Double.parseDouble(a) > 0.5 || Double.parseDouble(b) > 0.5){
        			if ( grafic.kindOfPhenotype == 5 ){
					dataTemp[i-1] = new Integer(i-1);
        			}
        			else
        			{
                		dataTemp[i-1] = new Integer(i);
        			}
       				k++;
			}

		}

		dataBun = new Object[k*nkat];
		k =0 ;
		for( int j =0 ; j < i; j++){
			if( dataTemp[j] != null){
				dataBun[k++] = dataTemp[j];
			}
		}
		if (nkat>1)
		{
			
			for( int j =k ; j < k*nkat; j++)
			{
				
				dataBun[j] = (Object)(new Integer(((Integer)dataBun[j-1]).intValue()+ 1));
				
			}
		}
	
		dataTemp = null;

		hypothesis = new Hipothesis(15,"First effect","Second effect","Hypothesis","hypothesis",dataBun);
		aditivity  = new Hipothesis(5,"First haplotype","Second haplotype","Deviation","deviation",dataBun);
		interaction = new Hipothesis(5,"Haplotype","Covariate","Interaction","interaction",2,dataBun,covar);
		homogeneity  = new Homogeneity(5,"Haplotype","1231","Hypothesis","hypothesis",dataBun);



		Border padding = BorderFactory.createEmptyBorder(20,20,5,20);

      
		jbtnApply  = new JButton(sbtnApply);
		jbtnCancel = new JButton(sbtnCancel);

		jbtnApply.setActionCommand(sbtnApply);
		jbtnApply.addActionListener(this);
		jbtnCancel.setActionCommand(sbtnCancel);
		jbtnCancel.addActionListener(this);


		JPanel butonsPanel = new JPanel();
		butonsPanel.add(jbtnApply,BorderLayout.CENTER);
		butonsPanel.add(jbtnCancel,BorderLayout.CENTER);

		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.addTab(sHypothesis, null, hypothesis, sHypothesisTip); 
		tabbedPane.addTab( sHomogeneity, null, homogeneity,sHomogeneityTip); 
		butonsPanel.setBackground(chair);
		butonsPanel.setForeground(chair);
		tabbedPane.setBackground(chair);
		jbtnApply.setBackground(saumon);
		jbtnCancel.setBackground(saumon);
		hypothesis.setBackground(chair);
		hypothesis.setForeground(chair);
		aditivity.setBackground(chair);
		interaction.setBackground(chair);
		homogeneity.setBackground(chair);
		
		 tabbedPane.addTab(sNonAdditivity, null,
	                       aditivity,
		                   sNonAdditivityTip); 
		 if ( covar.length != 0)
		 {
			 tabbedPane.addTab(sEnvironmentInteractions, null,
			                   interaction,
				               sEnvironmentInteractionsTip); 
		 }


		getContentPane().add(tabbedPane,BorderLayout.CENTER);	
		getContentPane().add(butonsPanel,BorderLayout.SOUTH);
		
		this.setSize(600,500);
	
			this.setBackground(saumon);

    }


    public void actionPerformed(ActionEvent e) {
		
       if ( sbtnCancel.equals(e.getActionCommand()) )
       {
		   setVisible(false);
       }
	   else if ( sbtnApply.equals(e.getActionCommand()) )
	   {
			saveDataToFile("paramData.thi");

			setCursor(hourglassCursor);
			grafic.t.run();
						StringBuffer sbuf = new StringBuffer();
						StringBuffer sBufText = new StringBuffer();
			String content;
			try
			{
				FileReader s;
				s = new FileReader("result.htm");
				BufferedReader b = new BufferedReader(s);
				content = b.readLine();
				
				while( content != null )
				{
					sbuf.append(content);
					content = b.readLine();
				}
				b.close();
				s.close();

			}
			catch ( IOException eIO	)
			{
				System.out.println(eIO.toString());
			}
			finally
			{}

				try
			{
				FileReader s;
				s = new FileReader("result.txt");
				BufferedReader b = new BufferedReader(s);
				content = b.readLine();
				
				while( content != null )
				{
					sBufText.append(content+"\r\n");
					content = b.readLine();
				}
				b.close();
				s.close();

			}
			catch ( IOException eIO	)
			{
				System.out.println(eIO.toString());
			}
			finally
			{}
			CustomDialog customDialog ;
			Object coaverCols2[] ;
			coaverCols2 = new Object[0];
			customDialog = new CustomDialog(datta,grafic,parent, "geisel",sbuf,new String(),coaverCols2,sBufText);

			setCursor(normalCursor);
			customDialog.setVisible(true);

	   }
    }


	
	public void saveDataToFile(String fileName)
	{
		String CR = "\r\n";
		StringBuffer sBuf = new StringBuffer();


		
		sBuf.append("y" + CR);	

		sBuf.append( hypothesis.saveData() );

		sBuf.append (homogeneity.saveData());

		sBuf.append( aditivity.saveData() );

		sBuf.append( interaction.saveData() );

		
		System.out.println(sBuf);
		FileWriter fwData;
		try
		{
			fwData = new FileWriter(fileName);
			fwData.write(sBuf.toString(), 0 , sBuf.toString().length());
			fwData.close();
		}
		catch (IOException e)
		{
			System.out.println(e.toString());
		}
	}

    public void propertyChange(PropertyChangeEvent e) {
		  
    }


    public void clearAndHide() {
        setVisible(false);
    }
}



	
	
 class SimpleDataModel extends AbstractTableModel {
        private String[] columnNames = {"Possible Effect"};
        private Object[][] data= {{new Integer(2)}};

		public SimpleDataModel(){
		}
		public SimpleDataModel(Object []left){
			data = new Object[left.length][1];
			for( int i =0 ; i < left.length; i++)
				data[i][0] = "Haplotype   " + ((Integer) left[i] ).toString();
		}

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return data.length;
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Object getValueAt(int row, int col) {
            return data[row][col];
        }
        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

        public boolean isCellEditable(int row, int col) {
                return false;
        }

        public void setValueAt(Object value, int row, int col) {
            data[row][col] = value;
            fireTableCellUpdated(row, col);
        }

    } 

	 class SimpleDataModel2D extends AbstractTableModel {
        private String[] columnNames = {"Possible Effect","Possible covariate"};
        private Object[][] data = {
            {new Integer(1),"V5"},
            {new Integer(2),"V6"},
            {new Integer(3),""},
            {new Integer(4),""},
        };
		public SimpleDataModel2D(){
		}
		public SimpleDataModel2D(Object []left,Object []right){
			int max;
			if ( left.length > right.length)
				max = left.length;
			else 
				max = right.length;
			data = new Object[max][2];
			for( int i =0 ; i< left.length; i++)
				data[i][0] =  "Haplotype  " + ((Integer) left[i] ).toString(); //left[i];
			for( int i =0 ; i< right.length; i++)
				data[i][1] = right[i];
		}
        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return data.length;
        }

        public String getColumnName(int col) {
            return columnNames[col];
        }

        public Object getValueAt(int row, int col) {
            return data[row][col];
        }
        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

        public boolean isCellEditable(int row, int col) {
                return false;
        }

        public void setValueAt(Object value, int row, int col) {
            data[row][col] = value;
            fireTableCellUpdated(row, col);
        }

    }

