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



class SimpleChooser extends JPanel  implements ActionListener 
	{
		JLabel lLeft,lRight,lTitle;
		JComboBox cLeft,cRight;
		public SimpleChooser(String l,String r,String t,Object []dataLeft,Object []dataRight){
			super ( new GridLayout(2,2) );
			lLeft = new JLabel(l);
			lRight = new JLabel(r);
			cLeft = new JComboBox(dataLeft);
			cRight = new JComboBox(dataRight);
			setBorder(BorderFactory.createTitledBorder(t) );
			add(lLeft);
			add(lRight);
			add(cLeft);
			add(cRight);
			cLeft.addActionListener(this);
			Color chair = new Color( 249, 230, 196 );
			Color saumon = new Color( 255, 204, 153 );
			lLeft.setBackground(chair);
			lRight.setBackground(chair);
			super.setBackground(chair);
			cLeft.setBackground(Color.white);
			cRight.setBackground(Color.white);
			
			
			
		}
		public SimpleChooser(String l,String r,String t,Object []data){
			super ( new GridLayout(2,2) );
			lLeft = new JLabel(l);
			lRight = new JLabel(r);
			cLeft = new JComboBox(data);
			cRight = new JComboBox(data);
			setBorder(BorderFactory.createTitledBorder(t) );
			add(lLeft);
			add(lRight);
			add(cLeft);
			add(cRight);
			cLeft.addActionListener(this);
			Color chair = new Color( 249, 230, 196 );
			Color saumon = new Color( 255, 204, 153 );
			cLeft.setBackground(Color.white);
			cRight.setBackground(Color.white);
			super.setBackground(chair);
		}
		
		
		public String getItems()
	{
        String s = (cRight.getSelectedItem()).toString();
        String t = (cRight.getSelectedItem()).toString();

        if( s.indexOf("(") != -1 ){
                t = s.substring(0,s.indexOf("(") );
        }
        return new String( ((Integer)cLeft.getSelectedItem()).toString() + " " + t );
} 

		
		
		public void actionPerformed(ActionEvent e) {

		}
	};

