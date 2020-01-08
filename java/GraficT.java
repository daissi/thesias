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
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.filechooser.*;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.table.*;

public class GraficT extends JPanel implements ActionListener, ItemListener  {
	final static boolean shouldFill = true;
	final static boolean shouldWeightX = true;
	final static boolean RIGHT_TO_LEFT = false;
	JFrame frame;
	JMenuBar menuBar;
	public JRadioButton v0Button, v1Button, v2Button, v3Button, v4Button, v5Button;
	ButtonGroup group;
	public JComboBox lpvCombo, idtimeCombo3_5, idtimeCombo4, numsxCombo, idoffsetCombo;
	JLabel l1,l2,l3;
	JButton runButton,advanceButton,exitButton;
	JCheckBox cbMissingGenotyp,cbPrintLDMatrix,cbEstimateHaplotypic,cbPrintR2,cbLink,cbOffset;
	JFileChooser fc;
	JTable table;
	MyTableModel tableModel;
	public boolean bMissingGenotyp, bPrintLDMatrix,bEstimateHaplotypic,bPrintR2,bLink,bOffset;
	public ListDemo lociChooser,covariables;
	public String rowstring,value;
	public Clipboard system;
	public StringSelection stsel;
	public int kindOfPhenotype;
	public boolean paramData;
	JScrollPane dataScrollPane; 
	public ParamTable paramTable;
	public thesiaslib t;
	private Cursor hourglassCursor;
	private Cursor normalCursor;
	private static Pattern pat1;
	private static Pattern pat2;
	private static Matcher mat1;
	private static Matcher mat2;
	public Color saumon;
	public Color chair; 
	public GraficT (  JFrame frame){
		super(new GridLayout(1,1));
		paramData = false;
		bMissingGenotyp = bPrintLDMatrix = bEstimateHaplotypic = bPrintR2 = bLink = bOffset =  false;
		this.frame = frame;
		addComponentsToPane();
		KeyStroke copy = KeyStroke.getKeyStroke(KeyEvent.VK_C,ActionEvent.CTRL_MASK,false);
		KeyStroke paste = KeyStroke.getKeyStroke(KeyEvent.VK_V,ActionEvent.CTRL_MASK,false);
		table.registerKeyboardAction(this,"Copy",copy,JComponent.WHEN_FOCUSED);
		table.registerKeyboardAction(this,"Paste",paste,JComponent.WHEN_FOCUSED);
		system = Toolkit.getDefaultToolkit().getSystemClipboard();
		hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
		normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
		fc = new JFileChooser("./");
		if (RIGHT_TO_LEFT) {setComponentOrientation(ComponentOrientation.RIGHT_TO_LEFT);}
		t = new thesiaslib();
		Color chair = new Color( 249, 230, 196 );
		frame.setBackground(chair);
		frame.setForeground(chair);			
	}
	public void addComponentsToPane(){
		Color chair = new Color( 249, 230, 196 );
		Color saumon = new Color( 255, 204, 153 );
		group = new ButtonGroup();
		v0Button = new JRadioButton("Null");v0Button.setActionCommand("v0Button");
		v0Button.setSelected(true);
		v0Button.setMnemonic(KeyEvent.VK_R);
		v1Button = new JRadioButton("Binary ");v1Button.setActionCommand("v1Button");
		v2Button = new JRadioButton("Quantitative");v2Button.setActionCommand("v2Button");
		v3Button = new JRadioButton("Survival");v3Button.setActionCommand("v3Button");
		v4Button = new JRadioButton("Matched C-C");	v4Button.setActionCommand("v4Button");
		v5Button = new JRadioButton("Categorical");	v5Button.setActionCommand("v5Button");
		group.add(v0Button);group.add(v1Button);	group.add(v2Button);
		group.add(v3Button);group.add(v4Button);group.add(v5Button);
		v0Button.addActionListener(this);v1Button.addActionListener(this);
		v2Button.addActionListener(this);v3Button.addActionListener(this);
		v4Button.addActionListener(this);v5Button.addActionListener(this); 
		JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
		splitPane.setResizeWeight(0.0);
		add(splitPane, BorderLayout.CENTER);
		splitPane.setBackground(chair);
		splitPane.setForeground(chair);
		v0Button.setBackground(chair);
		v1Button.setBackground(chair);
		v2Button.setBackground(chair);
		v3Button.setBackground(chair);
		v4Button.setBackground(chair);
		v5Button.setBackground(chair);
		JPanel tableLayout = new JPanel(new GridLayout(1,1) );
		tableLayout.setBorder(BorderFactory.createTitledBorder("Data table"));
		tableLayout.setBackground(chair);
		tableLayout.setForeground(chair);
		tableModel = new MyTableModel(7, 10);	
		table = new JTable(tableModel);
		JTableHeader tHeader;
		tHeader = table.getTableHeader();
		tHeader.setReorderingAllowed(false);
		tHeader.setResizingAllowed(true);
		tHeader.addMouseListener(new ColumnHeaderListener());
		table.setPreferredScrollableViewportSize(new Dimension(100, 120));
		table.setDragEnabled(true);
		table.setTableHeader(tHeader);
		tHeader.setBackground(saumon);
		table.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		table.setColumnSelectionAllowed(false);
		table.setRowSelectionAllowed(true);
		dataScrollPane = new JScrollPane(table);
		dataScrollPane.setRowHeaderView(computeRowHeader(7));
		tableLayout.add(dataScrollPane);
		splitPane.add(tableLayout);
		JPanel bottomLayout = new JPanel();
		GridBagLayout gridBagLayout = new GridBagLayout();
		GridBagConstraints cGridBag = new GridBagConstraints();
		bottomLayout.setLayout(gridBagLayout);
		bottomLayout.setBackground(chair);
		bottomLayout.setForeground(chair);
		JPanel checkPanel = new JPanel(new GridLayout(6, 0));
		checkPanel.setBorder(BorderFactory.createTitledBorder("Kind of phenotype") );
		checkPanel.add(v0Button); checkPanel.add(v1Button); checkPanel.add(v2Button);
		checkPanel.add(v3Button); checkPanel.add(v4Button); checkPanel.add(v5Button);	
		checkPanel.setBackground(saumon);
		JPanel paramVariables1 = new JPanel(new GridLayout(1,1));
		paramVariables1 .setBorder(BorderFactory.createTitledBorder("Selection of Variables") );
		paramVariables1.setBackground(saumon);
		paramVariables1.setForeground(saumon);
		JPanel paramVariables = new JPanel(new GridLayout(3, 2));
		paramVariables.setBackground(chair);
		paramVariables.setForeground(chair);
		JPanel temp = new JPanel(new GridLayout(6,1));
		temp.setBackground(saumon);
		temp.setForeground(saumon);
		cbMissingGenotyp = new JCheckBox("Deal with missing genotypic data");
		cbMissingGenotyp.addItemListener(this);
		cbPrintR2 =  new JCheckBox("Print R2 haplotypic 	information");
		cbPrintR2 .addItemListener(this);
		cbPrintLDMatrix  = new JCheckBox("Print LD Matrix");
		cbPrintLDMatrix.addItemListener(this);
		cbEstimateHaplotypic  = new JCheckBox("Estimation of haplotype effects");
		cbEstimateHaplotypic.addItemListener(this);
		temp.add(cbMissingGenotyp);
		temp.add(cbPrintR2);
		temp.add(cbPrintLDMatrix);
		temp.add(cbEstimateHaplotypic);
		JPanel tempBis1 = new JPanel(new GridLayout(1,2));
		JPanel tempBis2 = new JPanel(new GridLayout(1,2));
		temp.add(tempBis1);
		temp.add(tempBis2);
		cbLink  = new JCheckBox("X-Linked Analysis");
		cbLink.addItemListener(this);
		cbOffset = new JCheckBox("Offset");
		cbOffset.addItemListener(this);		
		tempBis1.add(cbLink);
		tempBis2.add(cbOffset);		
		cbMissingGenotyp.setBackground(chair);
		cbPrintR2.setBackground(chair);
		cbPrintLDMatrix.setBackground(chair);
		cbEstimateHaplotypic.setBackground(chair);
		cbLink.setBackground(chair);
		cbOffset.setBackground(chair);	
		String[] textForColomns = {"no value"};
		l1 = new JLabel("Phenotype variable");
		l1.setOpaque(true);
		l1.setBackground(chair); 
		lpvCombo = new JComboBox(textForColomns);		
		l2 = new JLabel("Time variable");
		l2.setOpaque(true);
		l2.setBackground(chair); 
		idtimeCombo3_5 = new JComboBox(textForColomns);
		l3 = new JLabel("Matching variable");
		l3.setOpaque(true);
		l3.setBackground(chair); 
		idtimeCombo4 = new JComboBox(textForColomns);
		numsxCombo = new JComboBox(textForColomns);
		idoffsetCombo = new JComboBox(textForColomns);
		lpvCombo.setBackground(Color.white);
		idtimeCombo3_5.setBackground(Color.white);
		idtimeCombo4.setBackground(Color.white);
		numsxCombo.setBackground(Color.white);
		idoffsetCombo.setBackground(Color.white);
		paramVariables.add(l1);
		paramVariables.add(lpvCombo);
		paramVariables.add(l2);
		paramVariables.add(idtimeCombo3_5);
		paramVariables.add(l3);
		paramVariables.add(idtimeCombo4);
		tempBis1.add(numsxCombo);
		tempBis2.add(idoffsetCombo);
		JPanel kindPanel = new JPanel(new GridLayout(2,1));
		kindPanel.setBackground(chair);
		kindPanel.add(checkPanel);
		kindPanel.add(paramVariables1);
		paramVariables1.add(paramVariables);
		cGridBag.gridx = 0;
		cGridBag.gridy = 0; 
		cGridBag.gridwidth  = 1;
		cGridBag.gridheight = 2;
		cGridBag.weightx =1.0;
		bottomLayout.add(kindPanel,cGridBag);
		cGridBag.gridx = 1;
		cGridBag.gridy = 0;
		cGridBag.gridwidth = 1;
		cGridBag.gridheight = 1;
		cGridBag.weightx =1.0;
		cGridBag.ipadx = 2;
		lociChooser = new ListDemo("Selection of Loci", this);	
		bottomLayout.add(lociChooser,cGridBag);
		lociChooser.setBackground(saumon);
		cGridBag.gridx = 2;
		cGridBag.gridy = 0;
		cGridBag.gridwidth = 1;
		cGridBag.gridheight = 1;
		cGridBag.weightx =1.0;
		covariables = new ListDemo("Selection of Covariates",null);
		covariables.setBackground(saumon);
		paramTable = new ParamTable();
		paramTable.setBackground(chair);
		cGridBag.gridx = 3;
		cGridBag.gridy = 0;
		cGridBag.gridwidth = 1;
		cGridBag.gridheight = 1;
		cGridBag.weightx =1.0;
		bottomLayout.add(paramTable,cGridBag);
		cGridBag.gridx = 1;
		cGridBag.gridy = 1;
		cGridBag.gridwidth = 2;
		cGridBag.gridheight = 1;
		bottomLayout.add(temp,cGridBag);
		cGridBag.gridx = 3;
		cGridBag.gridy = 1;
		cGridBag.gridwidth = 1;
		cGridBag.gridheight = 1;
		cGridBag.weightx = 0.001;
		cGridBag.weighty =.0;
		cGridBag.anchor = GridBagConstraints.SOUTH;
		JPanel butonPanel = new JPanel(new GridLayout(0,2));		
		runButton = new JButton("Run");
		butonPanel.add(runButton);
		runButton.setActionCommand("runButton");
		runButton.addActionListener(this);
		advanceButton = new JButton("Advance");
		advanceButton.setActionCommand("advanceButton");
		advanceButton.addActionListener(this);
		advanceButton.setVisible(false);
		exitButton = new JButton("Exit");
		exitButton.setActionCommand("exitButton");
		exitButton.addActionListener(this);
		butonPanel.add(exitButton);
		bottomLayout.add(butonPanel, cGridBag);
		exitButton.setBackground(saumon);
		runButton.setBackground(saumon);
		splitPane.add(bottomLayout);		
		setEditableText(0);
		v1Button.setEnabled(false);
		v2Button.setEnabled(false);
		v3Button.setEnabled(false);
		v4Button.setEnabled(false);
		v5Button.setEnabled(false);	
	}

	public JMenuBar addMenuBar(){
		menuBar = new JMenuBar();
		JMenu menu = new JMenu("Menu");
		JMenuItem menuItem = new JMenuItem("Open File");
		JMenuItem menuItem2 = new JMenuItem("About Thesias");
		Color saumon = new Color( 255, 204, 153 );
		menuBar.setBackground(saumon);
		menu.setBackground(saumon);
		menuItem.setBackground(saumon);
		menuItem2.setBackground(saumon);
		menuItem.setActionCommand("importMenu");
		menuItem.addActionListener(this);
		menuItem2.addActionListener(this);
		menu.add(menuItem);
		menu.add(menuItem2);
		menuBar.add(menu);
		return menuBar;
	}

	public void itemStateChanged(ItemEvent e){
		Object source = e.getItemSelectable();
		if (source == cbMissingGenotyp){
			bMissingGenotyp = ( bMissingGenotyp ) ? false : true;
		}
		else if ( source == cbPrintLDMatrix)
		{
			bPrintLDMatrix = ( bPrintLDMatrix ) ? false : true;
		}
		else if( source == cbEstimateHaplotypic )
		{
			bEstimateHaplotypic = ( bEstimateHaplotypic ) ? false : true;
			if( bEstimateHaplotypic )
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);

		}
		else if( source == cbPrintR2 ){
			bPrintR2 = ( bPrintR2 ) ? false : true;
		}

		else if(source == cbOffset){ 
			bOffset =(bOffset)? false : true;
			if(bOffset){idoffsetCombo.setEnabled(true);}
			else{idoffsetCombo.setEnabled(false);}		
		}
		
			else if(source == cbLink){ 
			bLink =(bLink)? false : true;
			if(bLink){numsxCombo.setEnabled(true);}
			else{numsxCombo.setEnabled(false);}		
		}
		

	}

	public void setEditableText(int i){
		kindOfPhenotype = i;
		switch( i )
		{
		case 0:
			lpvCombo.setEnabled(false);
			idtimeCombo3_5.setEnabled(false);
			idtimeCombo4.setEnabled(false);
			numsxCombo.setEnabled(false);
			idoffsetCombo.setEnabled(false);
			cbEstimateHaplotypic.setEnabled(false);
			cbLink.setEnabled(false);
			cbLink.setVisible(true);
			cbOffset.setEnabled(false);
			cbOffset.setVisible(true);
			paramTable.setVisible(true);
			paramTable.setParamTableEnabled(false);
			paramTable.reinitParamTable();
			lociChooser.setEnabled(false);
			covariables.setEnabled(false);
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbPrintR2.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
			lpvCombo.setSelectedIndex(0);
			idtimeCombo3_5.setSelectedIndex(0);
			idtimeCombo4.setSelectedIndex(0);
			numsxCombo.setSelectedIndex(0);
			idoffsetCombo.setSelectedIndex(0);		
			break;		
		case 1:
			lpvCombo.setEnabled(true);
			idtimeCombo3_5.setEnabled(false);
			idtimeCombo4.setEnabled(false);
			numsxCombo.setEnabled(false);
			idoffsetCombo.setEnabled(false);
			cbEstimateHaplotypic.setEnabled(true);
			cbLink.setVisible(true);
			cbLink.setEnabled(true);
			cbOffset.setVisible(true);
			cbOffset.setEnabled(true);
			paramTable.setVisible(true);
			paramTable.setParamTableEnabled(true);
			lociChooser.setEnabled(true);
			covariables.setEnabled(true);
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbPrintR2.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
			lpvCombo.setSelectedIndex(0);
			idtimeCombo3_5.setSelectedIndex(0);
			idtimeCombo4.setSelectedIndex(0);
			numsxCombo.setSelectedIndex(0);
			idoffsetCombo.setSelectedIndex(0);
			break;
		case 2:
			lpvCombo.setEnabled(true);
			idtimeCombo3_5.setEnabled(false);
			idtimeCombo4.setEnabled(false);
			numsxCombo.setEnabled(false);
			idoffsetCombo.setEnabled(false);
			cbEstimateHaplotypic.setEnabled(true);
			cbLink.setVisible(true);
			cbLink.setEnabled(true);
			cbOffset.setVisible(true);
			cbOffset.setEnabled(false);
			paramTable.setVisible(true);
			paramTable.setParamTableEnabled(true);
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbPrintR2.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
			lpvCombo.setSelectedIndex(0);
			idtimeCombo3_5.setSelectedIndex(0);
			idtimeCombo4.setSelectedIndex(0);
			numsxCombo.setSelectedIndex(0);
			idoffsetCombo.setSelectedIndex(0);
			break;
		case 3: 
			lpvCombo.setEnabled(true);
			idtimeCombo3_5.setEnabled(true);
			idtimeCombo4.setEnabled(false);
			numsxCombo.setEnabled(false);
			idoffsetCombo.setEnabled(false);
			cbEstimateHaplotypic.setVisible(true);
			cbEstimateHaplotypic.setEnabled(true);
			cbLink.setEnabled(false);
			cbLink.setVisible(true);
			cbOffset.setEnabled(false);
			cbOffset.setVisible(true);
			paramTable.setVisible(true);
			paramTable.setParamTableEnabled(true);
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbPrintR2.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
			lpvCombo.setSelectedIndex(0);
			idtimeCombo3_5.setSelectedIndex(0);
			idtimeCombo4.setSelectedIndex(0);
			numsxCombo.setSelectedIndex(0);
			idoffsetCombo.setSelectedIndex(0);
			break;
		case 4:
			lpvCombo.setEnabled(true);
			idtimeCombo3_5.setEnabled(false);
			idtimeCombo4.setEnabled(true);
			numsxCombo.setEnabled(false);
			idoffsetCombo.setEnabled(false);
			cbEstimateHaplotypic.setEnabled(true);
			cbLink.setEnabled(false);
			cbLink.setVisible(true);
			cbOffset.setEnabled(false);
			cbOffset.setVisible(true);
			paramTable.setVisible(true);
			paramTable.setParamTableEnabled(true);
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbPrintR2.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
			lpvCombo.setSelectedIndex(0);
			idtimeCombo3_5.setSelectedIndex(0);
			idtimeCombo4.setSelectedIndex(0);
			numsxCombo.setSelectedIndex(0);
			idoffsetCombo.setSelectedIndex(0);
			break;
		case 5:
			lpvCombo.setEnabled(true);
			idtimeCombo3_5.setEnabled(false);
			idtimeCombo4.setEnabled(false);
			numsxCombo.setEnabled(false);
			idoffsetCombo.setEnabled(false);
			cbEstimateHaplotypic.setEnabled(true);
			cbLink.setEnabled(false);
			cbLink.setVisible(true);
			cbOffset.setEnabled(false);
			cbOffset.setVisible(true);
			paramTable.setVisible(true);
			paramTable.setParamTableEnabled(true);
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbPrintR2.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
			lpvCombo.setSelectedIndex(0);
			idtimeCombo3_5.setSelectedIndex(0);
			idtimeCombo4.setSelectedIndex(0);
			numsxCombo.setSelectedIndex(0);
			idoffsetCombo.setSelectedIndex(0);
			break;
		}
	}
	
	public void actionPerformed(ActionEvent e) {
		if("About Thesias".equals(e.getActionCommand())){
			JOptionPane.showMessageDialog(this, "Testing Haplotype EffectS in Association Study (THESIAS)\nTHESIAS was developped by DA Tregouet,\nINSERM U525, Paris, France.\nFor any information, please contact DA Tregouet at\ndavid-alexandre.tregouet@inserm.fr", "About Thesias.",   JOptionPane.INFORMATION_MESSAGE,new ImageIcon("LogoThesias.jpg"));
}
		
		
		if( "v0Button".equals(e.getActionCommand()) ){
			setEditableText(0);
			advanceButton.setEnabled(false);
			lociChooser.setEnabled(true);
			lociChooser.ReInit(); 
			covariables.ReInit();
			covariables.setEnabled(false);
			paramTable.reinitParamTable();
			bMissingGenotyp = bPrintLDMatrix = bEstimateHaplotypic = bPrintR2 = bLink = bOffset =  false;
			cbMissingGenotyp.setSelected(false);
			cbPrintLDMatrix.setSelected(false);
			cbEstimateHaplotypic.setSelected(false);
			cbPrintR2.setSelected(false);
			cbLink.setSelected(false);
			cbOffset.setSelected(false);
		}
		
		else if( "v1Button".equals(e.getActionCommand()) ){
			setEditableText(1);
			covariables.setEnabled(true);
			if ( bEstimateHaplotypic ) 
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);
		}

		else if( "v2Button".equals(e.getActionCommand()) )
		{
			setEditableText(2);
			covariables.setEnabled(true);
			if ( bEstimateHaplotypic ) 
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);
		}
		else if( "v3Button".equals(e.getActionCommand()) ){
			setEditableText(3);
			covariables.setEnabled(true);
			if ( bEstimateHaplotypic ) 
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);
		}
		else if( "v4Button".equals(e.getActionCommand()) ){
			setEditableText(4);
			covariables.setEnabled(true);
			if ( bEstimateHaplotypic ) 
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);
		}

		else if( "v5Button".equals(e.getActionCommand()) ){
			setEditableText(5);
			covariables.setEnabled(true);
			if ( bEstimateHaplotypic ) 
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);

		}
		
		else if ( "importMenu".equals(e.getActionCommand()) ){
			int returnVal = fc.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {			
				File file = fc.getSelectedFile();
				lociChooser.setEnabled(true);			
				t = null;
				t = new thesiaslib();
				"v0Button".equals(e.getActionCommand());
				StringBuffer sbuf = new StringBuffer();
				String contentFile;
				int nbColumn=1;
				int nbLine=1;
				
				try{
					FileReader fImported = new FileReader(file);
					BufferedReader buffR1 = new BufferedReader(fImported);
					contentFile = buffR1.readLine();
					pat1 = Pattern.compile(" ");
					pat2 = Pattern.compile(";");
					mat1 = pat1.matcher(contentFile);
					mat2 = pat2.matcher(contentFile);					
					if (mat1.find()){contentFile = contentFile.replaceAll(" ", "\t");}
					if (mat2.find()){contentFile = contentFile.replaceAll(";", "\t");}
					StringTokenizer stTok = new StringTokenizer(contentFile,"\t");
					nbColumn = stTok.countTokens();
					int cpt;
					for (cpt = 0; buffR1.readLine() !=null;cpt++){
					}
					nbLine  =cpt+1;
					buffR1.close();
					fImported.close();
					fImported = null;
					buffR1 = null; 
					stTok = null;
				}
	
				catch (IOException excIO){
					System.out.println(excIO.toString());
				}
		
		tableModel = null;
		tableModel = new MyTableModel(nbLine, nbColumn);
		table.setModel(tableModel); 
		dataScrollPane.setRowHeaderView(computeRowHeader(nbLine));
		String content;
				try
				{
					int i,j,tj = 0; 
					Object d[][] = new Object[nbLine][nbColumn+1];
					FileReader s;
					s = new FileReader(file);
					BufferedReader b = new BufferedReader(s);
					content = b.readLine();
					pat1 = Pattern.compile(" ");
					pat2 = Pattern.compile(";");
					mat1 = pat1.matcher(content);
					mat2 = pat2.matcher(content);
					
					if (mat1.find()){content = content.replaceAll(" ", "\t");}
					if (mat2.find()){content = content.replaceAll(";", "\t");}
					StringTokenizer st;
					i = 0;
					while( content != null )
					{
					pat1 = Pattern.compile(" ");
					pat2 = Pattern.compile(";");
					mat1 = pat1.matcher(content);
					mat2 = pat2.matcher(content);
					
					if (mat1.find()){content = content.replaceAll(" ", "\t");}
					if (mat2.find()){content = content.replaceAll(";", "\t");}
					st = new StringTokenizer(content);
						j = 0;
						while (st.hasMoreTokens())
							d[i][j++] = st.nextToken();
						tj = j;
						i++;
						content = b.readLine();
					}
					b.close();
					s.close();

 b = null;
 s = null;
 st = null;
					tableModel.initData(d,i,tj);
					v1Button.setEnabled(false);
					v2Button.setEnabled(false);
					v3Button.setEnabled(false);
					v4Button.setEnabled(false);
					v5Button.setEnabled(false);
					lociChooser.ReInit();
					lociChooser.setTextForColumn(nbColumn);
					covariables.ReInit();
					paramTable.reinitParamTable();
					covariables.setTextForColumn(nbColumn);
					lpvCombo.setEnabled(false);
					idtimeCombo3_5.setEnabled(false);
					idtimeCombo4.setEnabled(false);
					numsxCombo.setEnabled(false);
					idoffsetCombo.setEnabled(false);
					setTextForColumn(nbColumn);
					paramTable.setVisible(true);
					v0Button.setSelected(true);					
					kindOfPhenotype = 0;
					cbMissingGenotyp.setSelected(false);
					cbPrintLDMatrix.setSelected(false);
					cbEstimateHaplotypic.setSelected(false);
					cbPrintR2.setSelected(false);
					cbLink.setSelected(false);
					cbOffset.setSelected(false);
					covariables.setEnabled(false);
				}
				catch ( IOException eIO	)
				{
					System.out.println(eIO.toString());
				}
			
			}
		}
		if ( e.getActionCommand().compareTo("exitButton") == 0 )
		{
			System.exit(0);
		}
		if (e.getActionCommand().compareTo("advanceButton") == 0)
		{}
		if( e.getActionCommand().compareTo("runButton")==0 )
		{
			if ( lociChooser.getsListSize() == 0)
			{
				return;
			}
			int i;
			int lpv = 0,ltv = 0,lmit = 0,lwv = 0, idoffset = 0, numx = 0;
			int idtime = 0;
			String paraFile = new String();
			Object data[][];

			if (v1Button.isSelected())
			{	
				if (!((String)lpvCombo.getSelectedItem()).equals(""))
				{
					lpv = Integer.parseInt( ((String)lpvCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
				
			}
			if (v2Button.isSelected())
			{
				if (!((String)lpvCombo.getSelectedItem()).equals(""))
				{
					lpv = Integer.parseInt( ((String)lpvCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for :Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
			}
			if (v3Button.isSelected()) 
			{
				if (!((String)lpvCombo.getSelectedItem()).equals(""))
				{
					lpv = Integer.parseInt( ((String)lpvCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}

				if (!((String)idtimeCombo3_5.getSelectedItem()).equals(""))
				{
					ltv = Integer.parseInt( ((String)idtimeCombo3_5.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
			}
			if (v4Button.isSelected())
			{
				if (!((String)lpvCombo.getSelectedItem()).equals(""))
				{
					lpv = Integer.parseInt( ((String)lpvCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
				if (!((String)idtimeCombo4.getSelectedItem()).equals(""))
				{
					lmit = Integer.parseInt( ((String)idtimeCombo4.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
			}

			if (v5Button.isSelected())
			{
				if (!((String)lpvCombo.getSelectedItem()).equals(""))
				{
					lpv = Integer.parseInt( ((String)lpvCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Need some value for Selection of variable !","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
			}

			
			if (cbLink.isEnabled() && cbLink.isSelected())
			{	
				if (!((String)numsxCombo.getSelectedItem()).equals(""))
				{
					numx = Integer.parseInt( ((String)numsxCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Please select a gender variable!","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
			}
			
			if (cbOffset.isEnabled() && cbOffset.isSelected())
			{	
				if (!((String)idoffsetCombo.getSelectedItem()).equals(""))
				{
					idoffset = Integer.parseInt( ((String)idoffsetCombo.getSelectedItem()).substring(1) );
				}
				else
				{
					JOptionPane.showMessageDialog(this, "Please select an offset variable!","WARNING", JOptionPane.WARNING_MESSAGE);
					return;
				}
			}
			
			
			else if( "v1Button".equals(e.getActionCommand()) ){
			setEditableText(1);
			covariables.setEnabled(true);
			if ( bEstimateHaplotypic ) 
				advanceButton.setEnabled(true);
			else 
				advanceButton.setEnabled(false);

		}

			if( kindOfPhenotype == 3 ) 
				idtime = ltv;
			else if (kindOfPhenotype == 4)
				idtime = lmit;

			if ( kindOfPhenotype > 0 )
			{
				paraFile = paramTable.serializeData();
			}


			Object objTemp[] = lociChooser.getElements();
			int lociCols[] = new int [lociChooser.getElements().length];
			for( int t = 0 ; t < lociCols.length ; t++ )
				lociCols[t] = Integer.parseInt( ((String)objTemp[t]).substring(1) );
			objTemp = covariables.getElements();
			int covarCols[];
			if (bLink)
			{
				covarCols = new int [covariables.getElements().length + 1];
				covarCols[0] = numx;
												
				for( int t = 1 ; t < covarCols.length ; t++ )
				{
					covarCols[t] = Integer.parseInt( ((String)objTemp[t-1]).substring(1) );
					
				}	
			}
			else
			{
				covarCols = new int [covariables.getElements().length];
				for( int t = 0 ; t < covarCols.length ; t++ )
				{
					covarCols[t] = Integer.parseInt( ((String)objTemp[t]).substring(1) );
				}
			}
			

			int col;

			data = tableModel.data;
			if (  ((String) data[0][0]).length() == 0  )
					return;
			v1Button.setEnabled(true);
			v2Button.setEnabled(true);
			v3Button.setEnabled(true);
			v4Button.setEnabled(true);
			v5Button.setEnabled(true);
			for(i = 0 ; i < data[0].length; i++)
				if (  ((String) data[0][i]).length() == 0  )
					break;
			col = i;
			for(i = 0; i < data.length ; i++)
			{	
				if( ((String) data[i][0]).length() == 0 )
					break;
			}
			data = new Object[i][col];
			for( i =0 ; i < data.length; i++)
				for(int j = 0 ; j < col ; j++)
					data[i][j] =  tableModel.data[i][j];
			
			writeDataToFile("datafile.txt", data, data.length, col);
			
			{
				FileWriter fwData;
				StringBuffer sBuf  = new StringBuffer();
				if( !bEstimateHaplotypic )
					sBuf.append("n" + "\r\n" + "n" + "\r\n" +  "n" + "\r\n" + "n" +"\r\n" +"n");
				else 
					sBuf.append("y" + "\r\n" + "n" + "\r\n" +  "n" + "\r\n" + "n" +"\r\n" +"n");
				try
				{
					fwData = new FileWriter("paramData.thi");
					fwData.write(sBuf.toString(), 0 , sBuf.toString().length());
					fwData.close();
				}
				catch (IOException ee)
				{
					System.out.println(ee.toString());
				}
			}

			dataThesias data3 = new dataThesias();
			data3.lociColsLength = lociCols.length;
			data3.bPrintLDMatrix = bPrintLDMatrix;
			data3.bPrintR2		 = bPrintR2;
			data3.bMissingGenotyp = bMissingGenotyp;
			data3.kindOfPhenotype = kindOfPhenotype;
			data3.lpv = lpv;
			data3.idtime= idtime;
			data3.lwv = idoffset;
			data3.covarColsLength = covarCols.length;
			data3.col = col;
			data3.initLoci(lociCols);
			data3.initCovar(covarCols);
			t.SetThesias(	"datafile.txt",
					lociCols.length,
					lociCols,
					(bPrintLDMatrix)  ? 1 : 0,
					(bMissingGenotyp) ? 1 : 0,
					(bPrintR2) ? 1 :0,
					kindOfPhenotype,
					lpv,
					idtime,
					(bOffset) ? 1 :0,
					idoffset,
					covarCols.length,
					covarCols,
					col,
					(bLink) ? 1 :0,
					numx);
			
			setCursor(hourglassCursor);
			t.run();
			System.gc();
			StringBuffer sbuf = new StringBuffer();
			StringBuffer sbufTest = new StringBuffer();
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
					sbufTest.append(content+"\r\n");
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
			coaverCols2 = new Object[covarCols.length];
			if ( covarCols.length > 0)
			{	
				for( int  k =0; k < covarCols.length ; k++ )
					coaverCols2[k] = (new Integer(covarCols[k])).toString();
			}
			if(bLink){
				
			
			System.out.println(covarCols.length);
				
			
			}
			
			
			
			customDialog = new CustomDialog(data3,this,frame, "geisel",sbuf,new String(),coaverCols2,sbufTest);
			if( !bEstimateHaplotypic )
				customDialog = new CustomDialog(data3,this,frame, "geisel",sbuf,new String(),coaverCols2,sbufTest);
			else 
				customDialog = new CustomDialog(data3,this,frame, "geisel",sbuf,paraFile,coaverCols2,sbufTest);

			setCursor(normalCursor);
			customDialog.setVisible(true);

			paramTable.importData();

		}
	    if (e.getActionCommand().compareTo("Copy")==0)
		{
         StringBuffer sbf=new StringBuffer();

         int numcols=table.getSelectedColumnCount();
         int numrows=table.getSelectedRowCount();
         int[] rowsselected=table.getSelectedRows();
         int[] colsselected=table.getSelectedColumns();
         if (!((numrows-1==rowsselected[rowsselected.length-1]-rowsselected[0] &&
                numrows==rowsselected.length) &&
(numcols-1==colsselected[colsselected.length-1]-colsselected[0] &&
                numcols==colsselected.length)))
         {
            JOptionPane.showMessageDialog(null, "Invalid Copy Selection",
                                          "Invalid Copy Selection",
                                          JOptionPane.ERROR_MESSAGE);
            return;
         }
         for (int i=0;i<numrows;i++)
         {
            for (int j=0;j<numcols;j++)
            {
sbf.append(table.getValueAt(rowsselected[i],colsselected[j]));
               if (j<numcols-1) sbf.append("\t");
            }
            sbf.append("\n");
         }
         stsel  = new StringSelection(sbf.toString());
         system = Toolkit.getDefaultToolkit().getSystemClipboard();
         system.setContents(stsel,stsel);
      }
      if (e.getActionCommand().compareTo("Paste")==0)
      {

          int startRow=(table.getSelectedRows())[0];
          int startCol=(table.getSelectedColumns())[0];
          try
          {
             String trstring= (String)(system.getContents(this).getTransferData(DataFlavor.stringFlavor));

             StringTokenizer st1=new StringTokenizer(trstring,"\n");
             for(int i=0;st1.hasMoreTokens();i++)
             {
                rowstring=st1.nextToken();
                StringTokenizer st2=new StringTokenizer(rowstring,"\t");
                for(int j=0;st2.hasMoreTokens();j++)
                {
                   value=(String)st2.nextToken();
                   if (startRow+i< table.getRowCount()  &&
                       startCol+j< table.getColumnCount())
                       table.setValueAt(value,startRow+i,startCol+j);

               }
            }
         }
         catch(Exception ex){ex.printStackTrace();}
      }
	}

	public void writeDataToFile(String name, Object [][]data, int d1, int d2)
	{
		FileWriter fwData;
		try
		{
			fwData = new FileWriter(name);
			StringBuffer buf;
			for( int i = 0 ; i < d1; i++){
				buf = new StringBuffer();
				for( int j = 0; j < d2-1; j++){
					buf.append( data[i][j] + "\t" );
				}
				buf.append( data[i][d2-1]);
				buf.append("\r\n");
				fwData.write(buf.toString(), 0 , buf.toString().length());

			}
			fwData.close();
			

		}
		catch (IOException e)
		{
			System.out.println(e.toString());
		}
	}
	

	public void setTextForColumn(int nbColumn){
		lpvCombo.removeAllItems();
		idtimeCombo3_5.removeAllItems();
		idtimeCombo4.removeAllItems();
		
		numsxCombo.removeAllItems();
		idoffsetCombo.removeAllItems();
			
		for (int j=0;j<1;j++){	
			String valCols = "";
			lpvCombo.addItem(valCols);
			idtimeCombo3_5.addItem(valCols);
			idtimeCombo4.addItem(valCols);
			numsxCombo.addItem(valCols);
			idoffsetCombo.addItem(valCols);
			j=1;
			
			for (int i=j;i<=nbColumn;i++){
				valCols = "V"+i;
				lpvCombo.addItem(valCols);
				idtimeCombo3_5.addItem(valCols);
				idtimeCombo4.addItem(valCols);
				numsxCombo.addItem(valCols);
				idoffsetCombo.addItem(valCols);
			}
		}
		
		String rien = " ";
		
	}


	public JList computeRowHeader(int nbLine)
	{	
		String[] headerTemp = new String[nbLine];		
		
		for (int i=0; i<nbLine; i++)  
		{

			headerTemp[i] = "" + (i+1);
		}
		final String[] finalHeaders = headerTemp;
		headerTemp = null;
		ListModel lm = new AbstractListModel() 
		{
			String[] headers = finalHeaders;
			public int getSize()
			{ 
				return headers.length;
			}
			public Object getElementAt(int index) 
			{
				return headers[index];
			}
		};
		
		JList rowHeader = new JList(lm);
		rowHeader.setFixedCellWidth(50);
		rowHeader.setFixedCellHeight(table.getRowHeight());
		rowHeader.setCellRenderer(new RowHeaderRenderer(table));

		return rowHeader;
	}

    public static void createAndShowGUI() {

JFrame.setDefaultLookAndFeelDecorated(false);

        JFrame frame = new JFrame("Thesias");
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		GraficT d = new GraficT(frame);
		frame.setJMenuBar(d.addMenuBar());
		d.setOpaque(true); 

		frame.setContentPane(d);

        frame.pack();

	frame.setSize(1024,600);
        frame.setVisible(true);
    }

    public static void main(String[] args) {
       javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {	
                createAndShowGUI();
            }
        });
    }
}

class MyTableModel extends AbstractTableModel {
	int rows, cols;
        private String[] columnNames;
        public Object[][] data;
	public MyTableModel( int row, int col){
		rows = row;
		cols = col;
		int i, j ;
		data = new Object[rows][cols];
		for( i =0 ; i < rows ; i++)
			for (j = 0 ; j < cols ; j++ ){
				data[i][j] = new String();
			}
			columnNames = new String[cols];
			for( i = 1; i <= cols ; i++)
				columnNames[i-1] = new String("V" + i);
	}

	public void resetData(){
		for( int i =0 ; i < rows ; i++)
			for (int j  = 0 ; j < cols ; j++ ){
				data[i][j] = new String();
			}
	}
	
	public void initData(Object d[][], int a, int b){
		for( int i = 0 ; i < a ; i++)
			for (int j  = 0 ; j < b ; j++ ){
				data[i][j] = (String)d[i][j];
			}
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
		if( data == null)
			return null;
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

class ListDemo extends JPanel implements ListSelectionListener {
	
	private JList list;
	public DefaultListModel listModel;
	private GraficT gT;
	
	private static final String hireString = "Add";
	private static final String fireString = "Remove";
	private JButton fireButton, hireButton;
	private JTextField employeeName;
	
	private JComboBox employeeNameCombo;
	private boolean multipleValues ;

	private String[] textForColomns = {""};
	
	
	
	public void setEnabled(boolean d){
		fireButton.setEnabled(d);	
		employeeName.setEnabled(d);
		employeeNameCombo.setEnabled(d);
		list.setEnabled(d);
		employeeNameCombo.setBackground(Color.white);
		
	}

	public ListDemo(String title, GraficT currentGT) {
		super(new BorderLayout());
		gT = currentGT;
		employeeNameCombo = new  JComboBox(textForColomns);
		multipleValues = false;
		listModel = new DefaultListModel();
		list = new JList(listModel);
		list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		list.addListSelectionListener(this);
		list.setVisibleRowCount(5);
		JScrollPane listScrollPane = new JScrollPane(list);
		hireButton = new JButton(hireString);
		HireListener hireListener = new HireListener(hireButton);
		hireButton.setActionCommand(hireString);
		hireButton.addActionListener(hireListener);
		hireButton.setEnabled(true);
		fireButton = new JButton(fireString);
		fireButton.setActionCommand(fireString);
		fireButton.addActionListener(new FireListener());
		fireButton.setEnabled(false);
		employeeName = new JTextField(5);
		employeeName.addActionListener(hireListener);
		employeeName.getDocument().addDocumentListener(hireListener);
		employeeNameCombo.addActionListener(hireListener);
		JPanel buttonPane = new JPanel();
		buttonPane.setLayout(new BoxLayout(buttonPane, BoxLayout.LINE_AXIS));
		buttonPane.add(fireButton);
		buttonPane.add(Box.createHorizontalStrut(5));
		buttonPane.add(Box.createHorizontalStrut(5));
		buttonPane.add(employeeNameCombo);
		buttonPane.add(hireButton);
		buttonPane.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		Color saumon = new Color( 255, 204, 153 );
		Color chair = new Color( 249, 230, 196 );
		employeeNameCombo.setBackground(Color.white);
		hireButton.setBackground(saumon);
		fireButton.setBackground(saumon);
		buttonPane.setBackground(chair);
		JLabel listTitle = new JLabel(title);
		add(listTitle, BorderLayout.NORTH);
		add(listScrollPane, BorderLayout.CENTER);
		add(buttonPane, BorderLayout.PAGE_END);
	}
	public void setTextForColumn(int nbColumn){
		employeeNameCombo.removeAllItems();
		for (int i=1;i<=nbColumn;i++){
			String valCols = "V"+i;
			employeeNameCombo.addItem(valCols);
		}
	}
	
	public void ReInit(){
		int index = list.getSelectedIndex();
		int size = listModel.getSize();
		size--;
		while( size != -1  ){
			listModel.remove(size--);	
		}
		fireButton.setEnabled(false); 
	}
	
	public int getsListSize(){
		return listModel.getSize();
	}


	public ListDemo(boolean multipleVal) {
        super(new BorderLayout());
		employeeNameCombo = new  JComboBox(textForColomns);
		multipleValues = multipleVal;
		listModel = new DefaultListModel();
		list = new JList(listModel);
		list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		list.addListSelectionListener(this);
		list.setVisibleRowCount(5);
		JScrollPane listScrollPane = new JScrollPane(list);
		JButton hireButton = new JButton(hireString);
		HireListener hireListener = new HireListener(hireButton,multipleVal);
		hireButton.setActionCommand(hireString);
		hireButton.addActionListener(hireListener);
		hireButton.setEnabled(true);
		fireButton = new JButton(fireString);
		fireButton.setActionCommand(fireString);
		fireButton.addActionListener(new FireListener());
		fireButton.setEnabled(false);
		employeeName = new JTextField(5);
		employeeName.addActionListener(hireListener);
		employeeName.getDocument().addDocumentListener(hireListener);
		employeeNameCombo.addActionListener(hireListener);
		JPanel buttonPane = new JPanel();
		buttonPane.setLayout(new BoxLayout(buttonPane, BoxLayout.LINE_AXIS));
		buttonPane.add(fireButton);
		buttonPane.add(Box.createHorizontalStrut(5));
		buttonPane.add(new JSeparator(SwingConstants.VERTICAL));
		buttonPane.add(Box.createHorizontalStrut(5));
		buttonPane.add(employeeNameCombo);
		buttonPane.add(hireButton);
		buttonPane.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		add(listScrollPane, BorderLayout.CENTER);
		add(buttonPane, BorderLayout.PAGE_END);
	}
	
	public Object[] getElements(){
		int i = listModel.getSize();
		Object data[] = new Object[i];
		for( int j =0 ; j < i ; j ++)
			data[j] = listModel.getElementAt(j);
		return data;
	 }

	class HireListener implements ActionListener, DocumentListener {
		private boolean alreadyEnabled = false;
		private JButton button;
		private boolean multipleValues;

		public HireListener(JButton button) {
			this.button = button;
			multipleValues = false;
		}
		public HireListener(JButton button,boolean multipleVal) {
		this.button = button;
		multipleValues = multipleVal;
        	}


        public void actionPerformed(ActionEvent e) {

			if( "Add".equals(e.getActionCommand()) &&
			    employeeNameCombo.isEnabled()){
	            String name = (String) employeeNameCombo.getSelectedItem() ;

				if (name.equals("") || alreadyInList(name)) {
					Toolkit.getDefaultToolkit().beep();
					employeeNameCombo.requestFocusInWindow();

					return;
				}

				int index = list.getSelectedIndex(); 
				if (index == -1) { 
					index = 0;
				} else {           
					index++;
				}
				listModel.insertElementAt(name, index);

				employeeName.requestFocusInWindow();

				list.setSelectedIndex(index);
				list.ensureIndexIsVisible(index);
				
				if (gT != null)
                                {
					gT.v0Button.setSelected(true);
					gT.lpvCombo.setEnabled(false);
					gT.idtimeCombo3_5.setEnabled(false);
					gT.idtimeCombo4.setEnabled(false);
					gT.numsxCombo.setEnabled(false);
					gT.idoffsetCombo.setEnabled(false);
					gT.cbEstimateHaplotypic.setEnabled(false);
					gT.cbLink.setEnabled(false);
					gT.cbOffset.setEnabled(false);
					gT.v1Button.setEnabled(false);
					gT.v2Button.setEnabled(false);
					gT.v3Button.setEnabled(false);
					gT.v4Button.setEnabled(false);
					gT.v5Button.setEnabled(false);
					gT.kindOfPhenotype = 0;
					gT.paramTable.reinitParamTable();
					gT.paramTable.setParamTableEnabled(false);
					gT.covariables.ReInit();
					gT.covariables.setEnabled(false);
				}
			
			}
			
        }

        protected boolean alreadyInList(String name) {
			if( multipleValues )
				return false;
            return listModel.contains(name);
        }


        public void insertUpdate(DocumentEvent e) {
            enableButton();
        }

 
        public void removeUpdate(DocumentEvent e) {
            handleEmptyTextField(e);
        }

        public void changedUpdate(DocumentEvent e) {
            if (!handleEmptyTextField(e)) {
                enableButton();
            }
        }

        private void enableButton() {
            if (!alreadyEnabled) {
                button.setEnabled(true);
            }
        }

        private boolean handleEmptyTextField(DocumentEvent e) {
            if (e.getDocument().getLength() <= 0) {
                button.setEnabled(false);
                alreadyEnabled = false;
                return true;
            }

			return false;
        }
    }

    class FireListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {

            int index = list.getSelectedIndex();
            listModel.remove(index);

            int size = listModel.getSize();

            if (size == 0) { 
                fireButton.setEnabled(false);
            } else { 
                if (index == listModel.getSize()) {
                    
                    index--;
                }

                list.setSelectedIndex(index);
                list.ensureIndexIsVisible(index);
            }

		if (gT != null){
			gT.v0Button.setSelected(true);
			gT.lpvCombo.setEnabled(false);
			gT.idtimeCombo3_5.setEnabled(false);
			gT.idtimeCombo4.setEnabled(false);
			gT.numsxCombo.setEnabled(false);
			gT.idoffsetCombo.setEnabled(false);
			gT.cbEstimateHaplotypic.setEnabled(false);
			gT.cbLink.setEnabled(false);
			gT.cbOffset.setEnabled(false);
			gT.v1Button.setEnabled(false);
			gT.v2Button.setEnabled(false);
			gT.v3Button.setEnabled(false);
			gT.v4Button.setEnabled(false);
			gT.v5Button.setEnabled(false);
			gT.kindOfPhenotype = 0;
			gT.paramTable.reinitParamTable();
			gT.paramTable.setParamTableEnabled(false);
			gT.covariables.ReInit();
			gT.covariables.setEnabled(false);
		}	
        }
    }



    public void valueChanged(ListSelectionEvent e) {
	if (e.getValueIsAdjusting() == false) {

            if (list.getSelectedIndex() == -1) {

                fireButton.setEnabled(false);

            } else {

                fireButton.setEnabled(true);
            }
        }
    }

    private static void createAndShowGUI() {

	JFrame.setDefaultLookAndFeelDecorated(false);

        JFrame frame = new JFrame("ListDemo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	JComponent newContentPane = new ListDemo("",null);
	
	newContentPane.setOpaque(true); 
        frame.setContentPane(newContentPane);


        frame.pack();
        frame.setVisible(true);
    }

}

class RowHeaderRenderer extends JLabel implements ListCellRenderer
{
   RowHeaderRenderer(JTable table) 
   {
	JTableHeader header = table.getTableHeader();
	setOpaque(true);
	setBorder(UIManager.getBorder("TableHeader.cellBorder"));
	setHorizontalAlignment(CENTER);
	setForeground(header.getForeground());
	setBackground(header.getBackground());
	setFont(header.getFont());
   }

   public Component getListCellRendererComponent(JList list,
						 Object value, 
						 int index, 
						 boolean isSelected, 
						 boolean cellHasFocus) 
   {
	setText((value == null) ? "" : value.toString());
	return this;
   }
}


class dataThesias{
	public int lociColsLength;
	public int lociCols[];
	public boolean bPrintLDMatrix;
	public boolean bPrintR2;
	public boolean bMissingGenotyp;
	public int kindOfPhenotype;
	public int lpv;
	public int idtime;
	public int lwv;
	public int covarColsLength;
	public int covarCols[];
	public int col;

	
	
	
	public void initLoci(int []loc){
		lociCols = new int [loc.length];
			for( int t = 0 ; t < loc.length ; t++ )
				lociCols[t] = loc[t];
	}
	
	public void initCovar(int []loc){
		covarCols = new int [loc.length];
		for( int t = 0 ; t < loc.length ; t++ )
			covarCols[t] = loc[t];
	}	
};


