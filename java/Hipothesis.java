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

	class Hipothesis extends JPanel  implements ActionListener 
	{
		String left,right;
		String title;
		SimpleChooser sch[];
		JTable tabel;
		AbstractTableModel dataModel;
		int nrRows;
		JButton addButton,removeButton;
		int visibleHypothesis;
		int tip;

		public Hipothesis(int maxRows,String firstEff,String secondEff,String mainTitle,String sButton,Object []left){
			super(new GridBagLayout());
			this.tip = 1;
			nrRows = maxRows;
			dataModel = new SimpleDataModel(left);
			makeAndDisplay(firstEff,secondEff,mainTitle,sButton,left,left);
		}

		public Hipothesis(int maxRows,String firstEff,String secondEff,String mainTitle,String sButton,int tip,Object []left,Object []right){
			super(new GridBagLayout());
			this.tip = 2;
			nrRows = maxRows;
			Object n[] = new Object[right.length];
			for( int i = 0; i < right.length ; i++ ) {
				n[i] = new String( (i+1) + " " + "(V" + right[i].toString() + ")" );
			}

			dataModel = new SimpleDataModel2D(left,n);
			makeAndDisplay(firstEff,secondEff,mainTitle,sButton,left,n);
			
		}

		private void makeAndDisplay(String firstEff,String secondEff,String mainTitle,String sButton,Object []left,Object []right)
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
			sch = new SimpleChooser[nrRows];
			for( int i =0 ;  i < nrRows ; i++){
				if( tip == 2)
					sch[i] = new SimpleChooser(firstEff,secondEff,mainTitle + " " +(i+1),left,right);
				else
					sch[i] = new SimpleChooser(firstEff,secondEff,mainTitle + " " +(i+1),left);
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
			StringBuffer buf = new StringBuffer();
			if (visibleHypothesis == 0)
			{
				buf.append("n\r\n");
			}
			else
			{
				buf.append("y\r\n");
				buf.append(visibleHypothesis + "\r\n");
				for( int i = 0 ; i < visibleHypothesis ; i++ )
				{
					buf.append(sch[i].getItems() + "\r\n") ;
				}

			}

			return new String(buf);

		}

		
	};
