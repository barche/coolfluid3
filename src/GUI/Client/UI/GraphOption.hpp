// Copyright (C) 2010 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef GRAPHOPTION_HPP
#define GRAPHOPTION_HPP

#include <QTableWidget>
#include <QPushButton>
#include <QComboBox>
#include <QDebug>

#include "qwt/qwt_plot.h"

#include "GUI/Client/Core/NPlotXY.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace CF {
namespace GUI {
namespace ClientUI {

////////////////////////////////////////////////////////////////////////////////

/// @brief This class is used to set visuals options, generate functions and
///  set the curve to draw.
/// @author Wertz Gil
class GraphOption : public QWidget
{
    Q_OBJECT
public: //functions

    //!\ => multi Array
    /// constructor
    /// @param fcts Data of functions.
    /// @param fcts_label Name of functions.
    /// @param ptr_plot Pointer of the QwtPlot where we draw functions.
    /// @param parent The QWidget of this QWidget.
    GraphOption(ClientCore::NPlotXY::PlotDataPtr & fcts,
                std::vector<QString> & fcts_label,
                QwtPlot * ptr_plot,
                QWidget *parent = 0);

    //!\ => multi Array
    /// Set the data to show in the option tab.
    /// @param fcts Data of functions.
    /// @param fcts_label Name of functions.
    void set_data(ClientCore::NPlotXY::PlotDataPtr & fcts,std::vector<QString> & fcts_label);

    //!\ => multi Array
    /// Add a function in the function set with it name and formula.
    /// @param fct Data of the function.
    /// @param function_name Name of the function.
    /// @param formula Formula of the function.
    void add_data(std::vector<double> & fct, QString fct_label ,QString formula = "");


private: //functions

    /// Same that the slot.
    void  generate_function(QString name,QString fct);

private: //datas

  /// QwtPlot pointer, refer to the plot where we draw cures.
  QwtPlot * m_ptr_plot;

  //!\ => multi Array
  /// The curves data.
  ClientCore::NPlotXY::PlotDataPtr m_fcts;

  /// User function's formula  line input.
  QLineEdit * m_line_function;

  /// User function's name line input.
  QLineEdit * m_line_function_name;

  /// Line table.
  QTableWidget * m_line_table;

  /// Data table.
  QTableWidget * m_data_table;

  /// Generate function button.
  QPushButton * button_generate_function;

  /// Draw line button.
  QPushButton * button_draw;

  bool m_can_draw;

private slots: //functions - slots

    /// Function that user inserted QString to make a function.
    void generate_function();

    /// Add a line to draw.
    void add_line();

    /// Remove selected line(s).
    void remove_line();

    /// Set the check state to each line selected and redraw curve(s).
    /// @param check_state The checked state.
    void checked_changed(int check_state);

    /// Set the color to each line selected and redraw curve(s).
    /// @param color The selected color.
    /// @param raw The row where the widget is stored.
    void color_changed(QColor color,int raw);

    /// Set the curve style to each line selected and redraw curve(s).
    /// @param current_index Index of the curve style.
    /// @param raw The row where the widget is stored.
    void line_type_changed(int current_index);

    /// Draw curves according to the data and options, and set axis auto scale.
    void draw_and_resize();

    /// Draw curves according to the data and options.
    void draw_action();

    /// Clear line table selection, unselect all line.
    void clear_line_table_selection();

    /// Select all raw of line table.
    void select_all_line_table();

};

////////////////////////////////////////////////////////////////////////////////

} // ClientUI
} // GUI
} // CF

#endif // GRAPHOPTION_HPP