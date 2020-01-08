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

import java.awt.*;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.datatransfer.*;
import java.util.*;
import javax.swing.JDialog;
import java.io.*;


	class Homogeneity extends JPanel  implements ActionListener 
	{
		String left,right;
		String title;
		HomogenitySimpleChooser sch[];
		JTable tabel;
		SimpleDataModel dataModel;

		int nrRows;
		JButton addButton,removeButton;
		int visibleHypothesis;
		private int tip;

		public Homogeneity(int maxRows,String firstEff,String secondEff,String mainTitle,String sButton,Object []left){
			super(new GridBagLayout());
			tip = 1;
			nrRows = maxRows;
			dataModel = new SimpleDataModel(left);
			makeAndDisplay(firstEff,secondEff,mainTitle,sButton,left);

		}
		private void makeAndDisplay(String firstEff,String secondEff,String mainTitle,String sButton,Object []left)
		{
			visibleHypothesis = 0;
			GridBagConstraints cGridBag = new GridBagConstraints();
			cGridBag.gridx = 0;
			cGridBag.gridy = 0;
			cGridBag.gridwidth = 1;
			cGridBag.gridheight = 2;
			cGridBag.weightx = 0.2;
			cGridBag.weighty = 1;
			cGridBag.fill = GridBagConstraints.BOTH;
			
			tabel = new JTable(dataModel);
			

			JTableHeader tHeader;
			tHeader = tabel.getTableHeader();
			Color saumon = new Color( 255, 204, 153 );
			tHeader.setBackground(saumon);		

			
			JScrollPane scrollPane = new JScrollPane(tabel);
			add(scrollPane,cGridBag);
		   
			cGridBag.gridx = 1;
			cGridBag.gridy = 0;
			cGridBag.gridwidth = 1;
			cGridBag.gridheight = 2;
			cGridBag.weightx = 0.3;
			cGridBag.weighty = 1;

			JPanel panel = new JPanel(new GridLayout(nrRows,1) );
			sch = new HomogenitySimpleChooser[nrRows];

			for( int i =0 ;  i < nrRows ; i++){
				sch[i] = new HomogenitySimpleChooser(firstEff,secondEff,mainTitle + " " +(i+1),left);
				sch[i].setVisible(false);
				panel.add(sch[i]);
			}
			JScrollPane scrollPane2 = new JScrollPane ( panel) ;
			add(scrollPane2,cGridBag);	
			JPanel buttonsPan = new JPanel(new GridLayout(2,0) );

			addButton = new JButton("Add " + sButton );
			addButton.addActionListener(this);
			removeButton = new JButton("Remove " + sButton);
			removeButton.addActionListener(this);
			
			buttonsPan.add(addButton);
			buttonsPan.add(removeButton);

			Color chair = new Color( 249, 230, 196 );
			tabel.setBackground(chair);
			panel.setBackground(chair);
			addButton.setBackground(saumon);
			removeButton.setBackground(saumon);
			scrollPane2.setBackground(saumon);
			scrollPane.setBackground(saumon);
			buttonsPan.setBackground(chair);
			buttonsPan.setForeground(chair);
			
			
			cGridBag.gridx = 3;
			cGridBag.gridy = 0;
			cGridBag.gridwidth = 1;
			cGridBag.gridheight = 1;
			cGridBag.weightx = 0;
			cGridBag.weighty = 0;
			cGridBag.fill = GridBagConstraints.NONE;
			add(buttonsPan,cGridBag);
		}	

		public void actionPerformed(ActionEvent e) {
			if ( addButton ==  (JButton)e.getSource()  )
			{
				if ( visibleHypothesis == nrRows )
					return ;
				sch[visibleHypothesis].setVisible(true);
				visibleHypothesis++;
			}
			else if ( removeButton ==  (JButton)e.getSource() )
			{
				if ( visibleHypothesis-1 < 0 )
					return ;
				visibleHypothesis--;
				sch[visibleHypothesis].setVisible(false);
			}
		}
		public String saveData()
		{
			StringBuffer  buf = new StringBuffer();
			if(visibleHypothesis == 0)
			{
				buf.append("n" + "\r\n");
			}
			else
			{
					buf.append("y" + "\r\n" + visibleHypothesis + "\r\n");
					for( int k = 0 ; k < visibleHypothesis  ; k++ ){
						buf.append( sch[k].getNomber() + "\r\n");
					}
					for( int k = 0 ; k < visibleHypothesis  ; k++ ){
						buf.append( sch[k].getPerechi() );
					}

			}
			return new String(buf);
		}
		
	};

	class  HomogenitySimpleChooser extends JPanel  implements ActionListener 
	{
		JLabel lLeft,lRight,lTitle;
		JComboBox cLeft,cRight;
		VerySimpleChooser vsch[];
		static int maxPairs = 8;
		JButton addPair,removePair;
		int visiblePair = 0;
		public HomogenitySimpleChooser(String l,String r,String t,Object []data){
			GridBagLayout gridBagLayout = new GridBagLayout();
			GridBagConstraints cGridBag = new GridBagConstraints();
			this.setLayout(gridBagLayout);
			Color chair = new Color( 249, 230, 196 );
			Color saumon = new Color( 255, 204, 153 );
			this.setBackground(chair);
			visiblePair = 2;
			addPair = new JButton("Add");
			removePair = new JButton("Remove");
			addPair.addActionListener(this);
			removePair.addActionListener(this);
			JPanel butP = new JPanel(new GridLayout(1,2));
			JPanel vschP = new JPanel( new GridLayout(4,2));
			butP.add(addPair);
			butP.add(removePair);
			butP.setBackground(chair);
			vschP.setBackground(chair);
			addPair.setBackground(saumon);
			removePair.setBackground(saumon);
			vsch = new VerySimpleChooser[maxPairs];
			setBorder(BorderFactory.createTitledBorder(t) );
			vsch[0] = new VerySimpleChooser(1,l,r,t,data,data);
			vsch[1] = new VerySimpleChooser(2,l,r,t,data,data);
			vschP.add(vsch[0]);
			vschP.add(vsch[1]);
			for(int i = 2 ; i < maxPairs; i++){
					vsch[i] = new VerySimpleChooser(i+1,l,r,t,data,data);

					vsch[i].setVisible(false);
					vschP.add(vsch[i]);
			}
		cGridBag.gridx = 0;
		cGridBag.gridy = 0;
		cGridBag.gridwidth  = 2;
		cGridBag.gridheight = 1;	

		cGridBag.fill = GridBagConstraints.HORIZONTAL;
		cGridBag.weightx = 1;
		cGridBag.weighty = 0.0;

			add(butP,cGridBag);
		cGridBag.gridx = 0;
		cGridBag.gridy = 1;
		cGridBag.gridwidth  = 2;
		cGridBag.gridheight = 4;	

		cGridBag.weightx = 1;
		cGridBag.weighty = 0.0;
			add(vschP,cGridBag);


		}
		public void actionPerformed(ActionEvent e) {
			if ( addPair ==  (JButton)e.getSource()  )
			{
				if ( visiblePair == maxPairs )
					return ;
				vsch[visiblePair].setVisible(true);
				visiblePair++;
			}
			else if ( removePair ==  (JButton)e.getSource() )
			{
				if ( visiblePair-1 < 2 )
					return ;
				visiblePair--;
				vsch[visiblePair].setVisible(false);
			}
		}
		public String getNomber(){
			return (new Integer(visiblePair)).toString();
		}
		public String getPerechi()
		{
			StringBuffer buf = new StringBuffer();
			for(int i =0 ; i < visiblePair ; i++){
				buf.append(vsch[i].getData() + "\r\n");
			}
			return new String(buf);
		}
	};

	class VerySimpleChooser extends JPanel
	{
		JLabel lLeft,lRight,lTitle;
		JComboBox cLeft,cRight;
		public VerySimpleChooser(int i,String l,String r,String t,Object []dataLeft,Object []dataRight){
			super ( new GridLayout(1,2) );
			setBorder(BorderFactory.createTitledBorder("Pair " + i) );

			cLeft = new JComboBox(dataLeft);
			cRight = new JComboBox(dataRight);
			add(cLeft);
			add(cRight);
			Color chair = new Color( 249, 230, 196 );
			super.setBackground(chair);
			cLeft.setBackground(Color.white);
			cRight.setBackground(Color.white);
			
		}
		public String getData()
		{
			return new String( ((Integer)cLeft.getSelectedItem()).toString() + " " + ((Integer)cRight.getSelectedItem()).toString() );
		}
	}
