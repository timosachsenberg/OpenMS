// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Matthias Fuchs $
// --------------------------------------------------------------------------



#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>

#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/FeatureSQLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/SqliteConnector.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>

#include <OpenMS/SYSTEM/File.h>

#include <unordered_set>
#include <sqlite3.h>
#include <sstream>
#include <map>
#include <boost/algorithm/string/join.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <boost/algorithm/string.hpp>

#include <typeinfo>

// #include <Boost>

using namespace std;

namespace OpenMS
{ 
  namespace Sql = Internal::SqliteHelper;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                            helper functions                                                          //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // store functions
  // convert int of enum datatype to a struct with prefix and type 
  PrefixSQLTypePair enumToPrefix(const DataValue::DataType& dt)
  {
    // create label prefix according to type
    String prefix;
    String type;
    PrefixSQLTypePair pSTP;
    switch (dt)
    {
      case DataValue::STRING_VALUE:
        pSTP.prefix = "_S_";
        pSTP.sqltype = "TEXT";
        break;
      case DataValue::INT_VALUE:
        pSTP.prefix = "_I_";
        pSTP.sqltype = "INTEGER";
        break;
      case DataValue::DOUBLE_VALUE:
        pSTP.prefix = "_D_";
        pSTP.sqltype = "FLOAT";
        break;
      case DataValue::STRING_LIST:
        pSTP.prefix = "_SL_";
        pSTP.sqltype = "TEXT";
        break;
      case DataValue::INT_LIST:
        pSTP.prefix = "_IL_";
        pSTP.sqltype = "TEXT";
        break;
      case DataValue::DOUBLE_LIST:
        pSTP.prefix = "_DL_";
        pSTP.sqltype = "TEXT";
        break;
      case DataValue::EMPTY_VALUE:
        pSTP.prefix = "";
        pSTP.sqltype = "TEXT";
        break;
    }
    return pSTP;
  }

  String createTable(const String& table_name, const String& table_stmt)
  {
    String create_table_stmt = "CREATE TABLE " + table_name + " (" + table_stmt + ");";
    return create_table_stmt;
  }

  // get datatype
  DataValue::DataType getColumnDatatype(const  String& label)
  { 
    DataValue::DataType type;
    if (label.hasPrefix("_S_"))
    {
      type = DataValue::STRING_VALUE;
    } 
    else if (label.hasPrefix("_I_"))
    {
      type = DataValue::INT_VALUE;
    } 
    else if (label.hasPrefix("_D_"))
    {
      type = DataValue::DOUBLE_VALUE;      
    } 
    else if (label.hasPrefix("_SL_"))
    {
      type = DataValue::STRING_LIST;
    } 
    else if (label.hasPrefix("_IL_"))
    {
      type = DataValue::INT_LIST;
    } 
    else if (label.hasPrefix("_DL_"))
    {
      type = DataValue::DOUBLE_LIST;
    } else 
    {
      type = DataValue::EMPTY_VALUE;
      /*
      if (!label.hasPrefix("_"))
      {
        type = 6;
      }
      */
    }
    return type;
  }

  /*
    // get userparameter
    String getColumnName(const String& label)
    { 
      String label_userparam;
      if (label.hasPrefix("_S_"))
      {
        label_userparam = label.substr(3);
      } 
      else if (label.hasPrefix("_I_"))
      {
        label_userparam = label.substr(3);
      } 
      else if (label.hasPrefix("_D_"))
      {
        label_userparam = label.substr(3);
      } 
      else if (label.hasPrefix("_SL_"))
      {
        label_userparam = label.substr(4);
      } 
      else if (label.hasPrefix("_IL_"))
      {
        label_userparam = label.substr(4);
      } 
      else if (label.hasPrefix("_DL_"))
      {
        label_userparam = label.substr(4);
      } else 
      {
        label_userparam = label;
        
        //if (!label.hasPrefix("_"))
        //{
        //  datatype = 6;
        //}
        
      }
      return label_userparam;
    }
  */

  // write userparameter to feature
  void setUserParams(Feature& current_feature, String &column_name, const DataValue::DataType column_type, int i, sqlite3_stmt* stmt)
  {
    if (column_type == DataValue::STRING_VALUE)
    {
      column_name = column_name.substr(3);
      String value;
      Sql::extractValue<String>(&value, stmt, i);
      current_feature.setMetaValue(column_name, value); 
    } 
    else if (column_type == DataValue::INT_VALUE)
    {
      column_name = column_name.substr(3);
      int value = 0;
      Sql::extractValue<int>(&value, stmt, i);
      current_feature.setMetaValue(column_name, value); 
    } 
    else if (column_type == DataValue::DOUBLE_VALUE)
    {
      column_name = column_name.substr(3);
      double value = 0.0;
      Sql::extractValue<double>(&value, stmt, i);          
      current_feature.setMetaValue(column_name, value); 
    } 
    else if (column_type == DataValue::STRING_LIST)
    {
      column_name = column_name.substr(4);
      String value; 
      Sql::extractValue<String>(&value, stmt, i);

      StringList sl;
      // cut off "[" and "]""
      value = value.chop(1);
      value = value.substr(1);
      value.split(", ", sl);
      current_feature.setMetaValue(column_name, sl);
    } 
    else if (column_type == DataValue::INT_LIST)
    {
      column_name = column_name.substr(4);
      String value; //IntList value;
      Sql::extractValue<String>(&value, stmt, i); //IntList
      value = value.chop(1);
      value = value.substr(1);
      std::vector<String> value_list;
      IntList il = ListUtils::create<int>(value, ',');
      current_feature.setMetaValue(column_name, il);
    } 
    else if (column_type == DataValue::DOUBLE_LIST)
    {
      column_name = column_name.substr(4);
      String value; //DoubleList value;
      Sql::extractValue<String>(&value, stmt, i); //DoubleList
      value = value.chop(1);
      value = value.substr(1);          
      DoubleList dl = ListUtils::create<double>(value, ',');
      current_feature.setMetaValue(column_name, dl);
    } 
    else if (column_type == DataValue::EMPTY_VALUE)
    {
      String value;
      Sql::extractValue<String>(&value, stmt, i);
    }
  }

  // probe first feature in feature_map as template to set tables
  // check once for features, subordinates, dataprocessing, convexhull bboxes respectively 
  std::tuple<bool, bool, bool, bool, bool> getTables(const FeatureMap& feature_map)
  {
    bool features_switch, subordinates_switch, dataprocessing_switch, features_bbox_switch, subordinates_bbox_switch;
    features_switch = false;
    features_bbox_switch = false;
    subordinates_switch = false;
    subordinates_bbox_switch = false;
    dataprocessing_switch = false;
    if (feature_map.getDataProcessing().size() > 0)
    {
      dataprocessing_switch = true;
    }
    for (auto feature : feature_map)
    {
      features_switch = true;
      if (feature.getConvexHulls().size() > 0)
      {
        features_bbox_switch = true;
      }
      for (Feature sub : feature.getSubordinates())
      {
        subordinates_switch = true;
        if (sub.getConvexHulls().size() > 0)
        {
          subordinates_bbox_switch = true;
        }
        break;
      }
      break;
    }
    return std::make_tuple(features_switch, subordinates_switch, dataprocessing_switch, features_bbox_switch, subordinates_bbox_switch);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                store FeatureMap as SQL database                                                      //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // fitted snippet of FeatureXMLFile::load
  void FeatureSQLFile::write(const std::string& out_fm, const FeatureMap& feature_map) const
  {
    String path="/home/mkf/Development/OpenMS/src/tests/class_tests/openms/data/";
    String filename = path.append(out_fm);
    
    // delete file if present
    File::remove(filename);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                        variable declaration                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

    // declare table entries
    std::vector<String> feature_elements = {"ID", "RT", "MZ", "Intensity", "Charge", "Quality"};
    std::vector<String> feature_elements_types = {"INTEGER", "REAL", "REAL", "REAL", "INTEGER", "REAL"};
    
    std::vector<String> subordinate_elements = {"ID", "SUB_IDX" , "REF_ID", "RT", "MZ", "Intensity", "Charge", "Quality"};
    std::vector<String> subordinate_elements_types = {"INTEGER", "INTEGER", "INTEGER", "REAL", "REAL", "REAL", "INTEGER", "REAL"};

    std::vector<String> dataprocessing_elements = {"ID", "SOFTWARE", "SOFTWARE_VERSION", "DATA", "TIME", "ACTIONS"};
    std::vector<String> dataprocessing_elements_types = {"INTEGER" ,"TEXT" ,"TEXT" ,"TEXT" ,"TEXT" , "TEXT"};

    std::vector<String> feat_bounding_box_elements = {"REF_ID", "MIN_MZ", "MIN_RT", "MAX_MZ", "MAX_RT", "BB_IDX"};
    std::vector<String> feat_bounding_box_elements_types = {"INTEGER" ,"REAL" ,"REAL" ,"REAL" , "REAL", "INTEGER"};

    std::vector<String> sub_bounding_box_elements = {"ID", "REF_ID", "MIN_MZ", "MIN_RT", "MAX_MZ", "MAX_RT", "BB_IDX"};
    std::vector<String> sub_bounding_box_elements_types = {"INTEGER" ,"INTEGER" ,"REAL" ,"REAL" ,"REAL" , "REAL", "INTEGER"};

    // get used tables
    auto tables = getTables(feature_map);
    bool features_switch = std::get<0>(tables);
    bool subordinates_switch = std::get<1>(tables);
    bool dataprocessing_switch = std::get<2>(tables);
    bool features_bbox_switch = std::get<3>(tables);
    bool subordinates_bbox_switch = std::get<4>(tables);

    std::cout << features_switch << "\t" << subordinates_switch << "\t" <<  dataprocessing_switch << "\t" << features_bbox_switch << "\t" << subordinates_bbox_switch << std::endl; 

    // read feature_map userparameter values and store as key value map
    // pass all CommonMetaKeys by setting frequency to 0.0 to a set 
    std::set<String> common_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<FeatureMap, std::set<String> >(feature_map.begin(), feature_map.end(), 0.0);
    std::map<String, DataValue::DataType> map_key2type;

    for (auto feature : feature_map)
    {
      for (const String& key : common_keys)
      {
        if (feature.metaValueExists(key))
        {
          const DataValue::DataType& dt = feature.getMetaValue(key).valueType();
          map_key2type[key] = dt;
        }
      }
    }

    std::vector<String> dataproc_keys;
    const std::vector<DataProcessing> dataprocessing_userparams = feature_map.getDataProcessing();
    std::map<String, DataValue::DataType> dataproc_map_key2type;

    for (auto datap : dataprocessing_userparams)
    {
      datap.getKeys(dataproc_keys);
      for (auto key : dataproc_keys)
      {
        const DataValue::DataType& dt = datap.getMetaValue(key).valueType(); 
        dataproc_map_key2type[key] = dt;

        /*
        const DataValue::DataType& dt = datap.getMetaValue(key).valueType(); 
        if ((dt == DataValue::STRING_VALUE) || (dt == DataValue::STRING_LIST))
        {
          dataproc_map_key2type[key] = dataproc_map_key2type["\"" + key + "\""];
        }
        dataproc_map_key2type.erase(key);
        */

      }
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// build feature header for sql table                                                                                                  //
    /// build subordinate header as sql table                                                                                               //
    /// build dataprocessing header as sql table                                                                                            //
    /// build boundingbox header as sql table                                                                                               //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // prepare feature header, add (dynamic) part of user_parameter labels to header and header type vector 
    for (const String& key : common_keys)
    {
      // feature_elements vector with strings prefix (_TYPE_, _S_, _IL_, ...) and key  
      feature_elements.push_back(enumToPrefix(map_key2type[key]).prefix + key);
      // feature_elements vector with SQL TYPES 
      feature_elements_types.push_back(enumToPrefix(map_key2type[key]).sqltype);
    }

    // prepare subordinate header
    std::map<String, DataValue::DataType> subordinate_key2type;
    for (FeatureMap::const_iterator feature_it = feature_map.begin(); feature_it != feature_map.end(); ++feature_it)
    {
      for (std::vector<Feature>::const_iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
      { 
        std::vector<String> keys;
        sub_it->getKeys(keys);
        for (const String& key : keys)
        {
          subordinate_key2type[key] = sub_it->getMetaValue(key).valueType();
        }
      }
    }
    for (const auto& key2type : subordinate_key2type)
    {
      // subordinate_elements vector with strings prefix (_TYPE_, _S_, _IL_, ...) and key  
      subordinate_elements.push_back(enumToPrefix(subordinate_key2type[key2type.first]).prefix + key2type.first);
      // subordinate_elements_type vector with SQL TYPES 
      subordinate_elements_types.push_back(enumToPrefix(subordinate_key2type[key2type.second]).sqltype);
    }

    // prepare dataprocessing header
    for (const auto& key2type : dataproc_map_key2type)
    {
      // dataprocessing vector with strings prefix (_TYPE_, _S_, _IL_, ...) and key  
      dataprocessing_elements.push_back(enumToPrefix(dataproc_map_key2type[key2type.first]).prefix + key2type.first);
      // dataprocessing_elements_type vector with SQL TYPES 
      dataprocessing_elements_types.push_back(enumToPrefix(key2type.second).sqltype);
    }
    // test for illegal entries
    // catch :
    const std::vector<String> bad_sym = {"+", "_", "-", "?", "!", "*", "@", "%", "^", "&", "#", "=", "/", "\\", ":", "\"", "\'"};

    for (std::size_t idx = 0; idx != dataprocessing_elements.size(); ++idx)
    {
      // test for illegal SQL entries in dataprocessing_elements, if found replace
      for (const String& sym : bad_sym)
      {
        if (dataprocessing_elements[idx].hasSubstring(sym))
        // use String::quote
        {
          dataprocessing_elements[idx] = dataprocessing_elements[idx].quote();
          break;
        }
      }
    }

    // boundingbox is part of convexhull is part of feature and subordinates
    // points of MZ/RT-min-/maxima saved as quadrupel of doubles per convexhull number
    // preparation of boundingbox header atm needs no implementation





    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                     create database with empty tables                                                //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // construct SQL_labels for feature, subordinate, dataprocessing and boundingbox tables
    // create empty
    // 1. features table
    // vector sql_labels, concatenate element and respective SQL type
    std::vector<String> sql_labels_features = {};
    for (std::size_t idx = 0; idx != feature_elements.size(); ++idx)
    {
      sql_labels_features.push_back(feature_elements[idx] + " " + feature_elements_types[idx]);
    }
    // add PRIMARY KEY to first entry, assuming 1st column contains primary key
    sql_labels_features[0].append(" PRIMARY KEY");
    // add "NOT NULL" to all entries
    std::for_each(sql_labels_features.begin(), sql_labels_features.end(), [] (String &s) {s.append(" NOT NULL");});
    // create single string with delimiter ',' as TABLE header template
    String sql_stmt_features = ListUtils::concatenate(sql_labels_features, ",");


    // 2. subordinates table
    // vector sql_labels, concatenate element and respective SQL type
    std::vector<String> sql_labels_subordinates = {};
    for (std::size_t idx = 0; idx != subordinate_elements.size(); ++idx)
    {
      sql_labels_subordinates.push_back(subordinate_elements[idx] + " " + subordinate_elements_types[idx]);
    }
    // add PRIMARY KEY to first entry, assuming 1st column contains primary key
    sql_labels_subordinates[0].append(" PRIMARY KEY");
    // add "NOT NULL" to all entries
    std::for_each(sql_labels_subordinates.begin(), sql_labels_subordinates.end(), [] (String &s) {s.append(" NOT NULL");});
    // create single string with delimiter ',' as TABLE header template
    String sql_stmt_subordinates = ListUtils::concatenate(sql_labels_subordinates, ",");


    // 3. dataprocessing table
    std::vector<String> sql_labels_dataprocessing = {};
    for (std::size_t idx = 0; idx != dataprocessing_elements.size(); ++idx)
    {
      sql_labels_dataprocessing.push_back(dataprocessing_elements[idx] + " " + dataprocessing_elements_types[idx]);
    }
    // add PRIMARY KEY to first entry, assuming 1st column contains primary key
    sql_labels_dataprocessing[0].append(" PRIMARY KEY");
    // add "NOT NULL" to all entries
    std::for_each(sql_labels_dataprocessing.begin(), sql_labels_dataprocessing.end(), [] (String &s) {s.append(" NOT NULL");});
    // create single string with delimiter ',' as TABLE header template
    String sql_stmt_dataprocessing = ListUtils::concatenate(sql_labels_dataprocessing, ",");


    // 4. boundingbox table
    // 4.1 features
    std::vector<String> sql_labels_boundingbox = {};
    for (std::size_t idx = 0; idx != feat_bounding_box_elements.size(); ++idx)
    {
      sql_labels_boundingbox.push_back(feat_bounding_box_elements[idx] + " " + feat_bounding_box_elements_types[idx]);
    }
    // no PRIMARY KEY because of executeStatement's UNIQUE constraint for primary keys
    // add "NOT NULL" to all entries
    std::for_each(sql_labels_boundingbox.begin(), sql_labels_boundingbox.end(), [] (String &s) {s.append(" NOT NULL");});
    // create single string with delimiter ',' as TABLE header template
    String sql_stmt_feat_boundingbox = ListUtils::concatenate(sql_labels_boundingbox, ",");


    // 4. boundingbox table
    // 4.2 subordinates
    sql_labels_boundingbox = {};
    for (std::size_t idx = 0; idx != sub_bounding_box_elements.size(); ++idx)
    {
      sql_labels_boundingbox.push_back(sub_bounding_box_elements[idx] + " " + sub_bounding_box_elements_types[idx]);
    }
    // add PRIMARY KEY to first entry, assuming 1st column contains primary key
    //sql_labels_boundingbox[0].append(" PRIMARY KEY");
    // add "NOT NULL" to all entries
    std::for_each(sql_labels_boundingbox.begin(), sql_labels_boundingbox.end(), [] (String &s) {s.append(" NOT NULL");});
    // create single string with delimiter ',' as TABLE header template
    String sql_stmt_sub_boundingbox = ListUtils::concatenate(sql_labels_boundingbox, ",");


    // create empty SQL table stmt
    String features_table_stmt, subordinates_table_stmt, dataprocessing_table_stmt, feature_boundingbox_table_stmt, subordinate_boundingbox_table_stmt;
    // 1. features
    if (features_switch == true)
    {
      features_table_stmt = createTable("FEATURES_TABLE", sql_stmt_features);
    } 
    else if (features_switch == false)
    {
      features_table_stmt = "";
    }
    // 2. subordinates
    if (subordinates_switch == true)
    {
      subordinates_table_stmt = createTable("FEATURES_SUBORDINATES", sql_stmt_subordinates);
    } 
    else if (subordinates_switch == false)
    {
      subordinates_table_stmt = "";
    }
    // 3. dataprocessing
    if (dataprocessing_switch == true)
    {
      dataprocessing_table_stmt = createTable("FEATURES_DATAPROCESSING", sql_stmt_dataprocessing);
    } 
    else if (dataprocessing_switch == false)
    {
      dataprocessing_table_stmt = "";
    }
    // 4. boundingbox (features & subordinates)
    // features
    if (features_bbox_switch == true)
    {
      feature_boundingbox_table_stmt = createTable("FEATURES_TABLE_BOUNDINGBOX", sql_stmt_feat_boundingbox);
    } 
    else if (features_bbox_switch == false)
    {
      feature_boundingbox_table_stmt = "";
    }
    // subordinates
    if (subordinates_bbox_switch)
    {
      subordinate_boundingbox_table_stmt = createTable("SUBORDINATES_TABLE_BOUNDINGBOX", sql_stmt_sub_boundingbox);
    } 
    else if (subordinates_bbox_switch == false)
    {
      subordinate_boundingbox_table_stmt = "";
    }


    String create_sql = \
      features_table_stmt + \
      subordinates_table_stmt + \
      dataprocessing_table_stmt + \
      feature_boundingbox_table_stmt + \
      subordinate_boundingbox_table_stmt \
      ;    

    // catch SQL contingencies
    //create_sql = create_sql.substitute("'", "\:"); // SQL escape single quote in strings  
    // create SQL database, create empty SQL tables  see (https://github.com/OpenMS/OpenMS/blob/develop/src/openms/source/FORMAT/OSWFile.cpp)
    // Open connection to database
    SqliteConnector conn(filename);
    conn.executeStatement(create_sql);





































    // break into functions?
    // turn into line feed function for SQL database input step statement
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                     store FeatureMap data in tables                                                  //
    //                                                     1. features                                                                      //                                    
    //                                                     2. subordinates                                                                  //
    //                                                     3. dataprocessing                                                                //
    //                                                     4. boundingboxes features & subordinates                                         //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// 1. insert data of features table
    const String feature_elements_sql_stmt = ListUtils::concatenate(feature_elements, ","); 
    const String feat_bbox_elements_sql_stmt = ListUtils::concatenate(feat_bounding_box_elements, ","); 

    if (features_switch)
    {
      for (const Feature& feature : feature_map)
      {
        String line_stmt_features_table;
        std::vector<String> line;

        int64_t id = static_cast<int64_t>(feature.getUniqueId() & ~(1ULL << 63));
        line.push_back(id);
        line.push_back(feature.getRT());
        line.push_back(feature.getMZ());
        line.push_back(feature.getIntensity());
        line.push_back(feature.getCharge());
        line.push_back(feature.getOverallQuality());
        for (const String& key : common_keys)
        {
          String s = feature.getMetaValue(key);
          if (map_key2type[key] == DataValue::STRING_VALUE
          || map_key2type[key] == DataValue::STRING_LIST)
          {
            s = s.substitute("'", "''"); // SQL escape single quote in strings  
          }

          if (map_key2type[key] == DataValue::STRING_VALUE
          || map_key2type[key] == DataValue::STRING_LIST
          || map_key2type[key] == DataValue::INT_LIST
          || map_key2type[key] == DataValue::DOUBLE_LIST)
          {
            s = "'" + s + "'"; // quote around SQL strings (and list types)
          }
          line.push_back(s);
        }
        line_stmt_features_table =  "INSERT INTO FEATURES_TABLE (" + feature_elements_sql_stmt + ") VALUES (";
        line_stmt_features_table += ListUtils::concatenate(line, ",");
        line_stmt_features_table += ");";
        //store in features table
        conn.executeStatement(line_stmt_features_table);
      }
    }

    if (features_bbox_switch)
    {
      // fetch convexhull of feature for bounding box data
      for (const Feature& feature : feature_map)
      {
        /// insert boundingbox data to table FEATURES_TABLE_BOUNDINGBOX
        std::vector<String> feat_bbox_line;

        double min_MZ;
        double min_RT;
        double max_MZ;
        double max_RT;

        int bbox_idx = 0;
        int64_t id = static_cast<int64_t>(feature.getUniqueId() & ~(1ULL << 63));

        // add bbox entries of current convexhull
        for (Size b = 0; b < feature.getConvexHulls().size(); ++b)
        {
          //it->getConvexHulls()[i].expandToBoundingBox();

          //clear vector line for current bbox
          feat_bbox_line.clear();
          String line_stmt_features_table_boundingbox =  "INSERT INTO FEATURES_TABLE_BOUNDINGBOX (" + feat_bbox_elements_sql_stmt + ") VALUES (";
          bbox_idx = b;

          min_MZ = feature.getConvexHulls()[b].getBoundingBox().minX();
          min_RT = feature.getConvexHulls()[b].getBoundingBox().minY();
          max_MZ = feature.getConvexHulls()[b].getBoundingBox().maxX();
          max_RT = feature.getConvexHulls()[b].getBoundingBox().maxY();

          feat_bbox_line.push_back(id);
          feat_bbox_line.insert(feat_bbox_line.end(), {min_MZ,min_RT,max_MZ,max_RT});
          feat_bbox_line.push_back(bbox_idx);

          line_stmt_features_table_boundingbox += ListUtils::concatenate(feat_bbox_line, ",");
          //store String in FEATURES_TABLE_BOUNDINGBOX
          line_stmt_features_table_boundingbox += ");";
          conn.executeStatement(line_stmt_features_table_boundingbox);
        }
      }
    }

    if (subordinates_switch)
    {
      /// 2. insert data of subordinates table
      String line_stmt;
      const String subordinate_elements_sql_stmt = ListUtils::concatenate(subordinate_elements, ","); 
      //for (auto feature : feature_map)
      for (const Feature& feature : feature_map)
      {
        int sub_idx = 0;
        std::vector<String> line;
        for (const Feature& sub : feature.getSubordinates())
        {          
          line.push_back(static_cast<int64_t>(sub.getUniqueId() & ~(1ULL << 63)));
          // additional index value to preserve order of subordinates
          line.push_back(sub_idx);
          ++sub_idx;
          line.push_back(static_cast<int64_t>(feature.getUniqueId() & ~(1ULL << 63)));
          line.push_back(sub.getRT());
          line.push_back(sub.getMZ());
          line.push_back(sub.getIntensity());
          line.push_back(sub.getCharge());
          line.push_back(sub.getOverallQuality());
          for (const std::pair<String, DataValue::DataType> k2t : subordinate_key2type)
          {
            const String& key = k2t.first;
            if (sub.metaValueExists(key))
            {
              String s = sub.getMetaValue(k2t.first);
              if (k2t.second == DataValue::STRING_VALUE
              || k2t.second == DataValue::STRING_LIST)
              {
                s = s.substitute("'", "''"); // SQL escape single quote in strings  
              }

              if (k2t.second == DataValue::STRING_VALUE
              || k2t.second == DataValue::STRING_LIST
              || k2t.second == DataValue::INT_LIST
              || k2t.second == DataValue::DOUBLE_LIST)
              {
                s = "'" + s + "'"; // quote around SQL strings (and list types)
              }
              line.push_back(s);
            }
          }
          line_stmt =  "INSERT INTO FEATURES_SUBORDINATES (" + subordinate_elements_sql_stmt + ") VALUES (";
          line_stmt += ListUtils::concatenate(line, ",");
          line_stmt += ");";
          //store in subordinates table
          conn.executeStatement(line_stmt);
          //clear vector line
          line.clear();
        }
      }
    }

    if (subordinates_bbox_switch)
    {
      // fetch bounding box data of subordinate convexhull
      const String sub_bbox_elements_sql_stmt = ListUtils::concatenate(sub_bounding_box_elements, ","); 
      for (const Feature& feature : feature_map)
      {
        int64_t ref_id = static_cast<int64_t>(feature.getUniqueId() & ~(1ULL << 63));

        for (const Feature& sub : feature.getSubordinates())
        {
          int64_t id = static_cast<int64_t>(sub.getUniqueId() & ~(1ULL << 63));

          /// insert boundingbox data to table SUBORDINATES_TABLE_BOUNDINGBOX
          std::vector<String> bbox_line;

          double min_MZ;
          double min_RT;
          double max_MZ;
          double max_RT;

          int bbox_idx = 0;

          // add bbox entries of current convexhull
          for (Size b = 0; b < sub.getConvexHulls().size(); ++b)
          {
            //clear vector line for current bbox
            bbox_line.clear();
            String line_stmt_subordinates_table_boundingbox =  "INSERT INTO SUBORDINATES_TABLE_BOUNDINGBOX (" + sub_bbox_elements_sql_stmt + ") VALUES (";
            bbox_idx = b;

            min_MZ = sub.getConvexHulls()[b].getBoundingBox().minX();
            min_RT = sub.getConvexHulls()[b].getBoundingBox().minY();
            max_MZ = sub.getConvexHulls()[b].getBoundingBox().maxX();
            max_RT = sub.getConvexHulls()[b].getBoundingBox().maxY();

            bbox_line.push_back(id);
            bbox_line.push_back(ref_id);
            bbox_line.insert(bbox_line.end(), {min_MZ,min_RT,max_MZ,max_RT});
            bbox_line.push_back(bbox_idx);

            line_stmt_subordinates_table_boundingbox += ListUtils::concatenate(bbox_line, ",");
            line_stmt_subordinates_table_boundingbox += ");";
            conn.executeStatement(line_stmt_subordinates_table_boundingbox);
          }
        }
      }
    }

    if (dataprocessing_switch)
    {
      /// 3. insert data of dataprocessing table
      const String dataprocessing_elements_sql_stmt = ListUtils::concatenate(dataprocessing_elements, ","); 
      std::vector<String> dataproc_elems = {};
      std::vector<String> processing_action_names;
      

      std::cout << "################################################" << std::endl;
      std::cout << "\n" << static_cast<int64_t>(feature_map.getUniqueId() & ~(1ULL << 63)) << std::endl;
      std::cout << "\n" << feature_map.getUniqueId() << std::endl;
      std::cout << "################################################" << std::endl;
      
      dataproc_elems.push_back(static_cast<int64_t>(feature_map.getUniqueId() & ~(1ULL << 63)));
      //dataproc_elems.push_back(feature_map.getUniqueId());
      
      const std::vector<DataProcessing> dataprocessing = feature_map.getDataProcessing();
    
      for (auto datap : dataprocessing)
      {
        dataproc_elems.push_back(datap.getSoftware().getName());
        dataproc_elems.push_back(datap.getSoftware().getVersion());
        dataproc_elems.push_back(datap.getCompletionTime().getDate());
        dataproc_elems.push_back(datap.getCompletionTime().getTime());

        const std::set<DataProcessing::ProcessingAction>& proc_ac = datap.getProcessingActions();

        for (const DataProcessing::ProcessingAction& a : proc_ac)
        {
          const String& action = DataProcessing::NamesOfProcessingAction[int(a)];
          processing_action_names.push_back(action);
        }
        if (processing_action_names.size() >= 1)
        {
          dataproc_elems.push_back(ListUtils::concatenate(processing_action_names, ",")); 
        } else
        {
          dataproc_elems.push_back(processing_action_names[0]);
        }
      }
    
      // userparam entries
      // dynamic type resolution 
      for (auto datap : dataprocessing_userparams)
      {
        datap.getKeys(dataproc_keys);
        for (auto key : dataproc_keys)
        {
          dataproc_elems.push_back(datap.getMetaValue(key));
        }
      }

      // non-userparam entries
      // resolve SQL type from dataprocessing_elements_types
      for (std::size_t idx = 0; idx != dataprocessing_elements.size(); ++idx)
      {
        if (dataprocessing_elements_types[idx] == "TEXT")
        {
          dataproc_elems[idx] = "'" + dataproc_elems[idx] + "'"; 
        }
      }
      String line_stmt =  "INSERT INTO FEATURES_DATAPROCESSING (" + dataprocessing_elements_sql_stmt + ") VALUES (";
      line_stmt += ListUtils::concatenate(dataproc_elems, ",");
      line_stmt += ");";
      //store in dataprocessing table
      conn.executeStatement(line_stmt);
    }

  } // end of FeatureSQLFile::write


  // read functions
<<<<<<< HEAD

  ConvexHull2D readBBox_(sqlite3_stmt * stmt, int  column_nr)
=======
  ConvexHull2D readBBox_(sqlite3_stmt * stmt, int column_nr)
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
  {
    std::cout << "\n column_nr = " << column_nr << std::endl;

    double min_mz = 0.0;
    Sql::extractValue<double>(&min_mz, stmt, (column_nr + 1));
    double min_rt = 0.0;
    Sql::extractValue<double>(&min_rt, stmt, (column_nr + 2));
    double max_mz = 0.0;
    Sql::extractValue<double>(&max_mz, stmt, (column_nr + 3));
    double max_rt = 0.0;
    Sql::extractValue<double>(&max_rt, stmt, (column_nr + 4));

    ConvexHull2D hull;
    hull.addPoint({min_mz, min_rt});
    hull.addPoint({max_mz, max_rt});
    return hull;
  }

  // 
  Feature readSubordinate_(sqlite3_stmt * stmt, int column_nr, int cols_features, int cols_subordinates)
  {
    Feature subordinate;

    String id_string;
    Sql::extractValue<String>(&id_string, stmt, cols_features);
    std::istringstream iss(id_string);
    int64_t id_sub = 0;
    iss >> id_sub;
<<<<<<< HEAD

=======
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5

    double rt = 0.0;
    Sql::extractValue<double>(&rt, stmt, column_nr + 1);
    double mz = 0.0;
    Sql::extractValue<double>(&mz, stmt, column_nr + 2);
    double intensity = 0.0;
    Sql::extractValue<double>(&intensity, stmt, column_nr + 3);
    int charge = 0;
    Sql::extractValue<int>(&charge, stmt, column_nr + 4);
    double quality = 0.0;
    Sql::extractValue<double>(&quality, stmt, column_nr + 5);
<<<<<<< HEAD
    subordinate.setUniqueId(id_sub  );
=======
    subordinate.setUniqueId(id_sub);
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
    subordinate.setRT(rt);
    subordinate.setMZ(mz);
    subordinate.setIntensity(intensity);
    subordinate.setCharge(charge);
    subordinate.setOverallQuality(quality);

    // userparams
    column_nr = cols_features + 5;  // +5 = subordinate params
    for (int i = column_nr; i < cols_features + cols_subordinates  ; ++i) // subordinate userparams - 1 = index of last element in cols_subordinates
    {
      String column_name = sqlite3_column_name(stmt, i);
      int column_type = getColumnDatatype(column_name);


      if (column_type == DataValue::STRING_VALUE)
      {
        column_name = column_name.substr(3);
        String value;
        Sql::extractValue<String>(&value, stmt, i);
        subordinate.setMetaValue(column_name, value); 
        continue;
      } 
      else if (column_type == DataValue::INT_VALUE)
      {
        column_name = column_name.substr(3);
        int value = 0;
        Sql::extractValue<int>(&value, stmt, i);
        subordinate.setMetaValue(column_name, value); 
        continue;
      } 
      else if (column_type == DataValue::DOUBLE_VALUE)
      {
        column_name = column_name.substr(3);
        double value = 0.0;
        Sql::extractValue<double>(&value, stmt, i);          
        subordinate.setMetaValue(column_name, value); 
        continue;
      } 
      else if (column_type == DataValue::STRING_LIST)
      {
        column_name = column_name.substr(4);
        String value; 
        Sql::extractValue<String>(&value, stmt, i);

        StringList sl;
        // cut off "[" and "]""
        value = value.chop(1);
        value = value.substr(1);
        value.split(", ", sl);
        subordinate.setMetaValue(column_name, sl);
        continue;
      } 
      else if (column_type == DataValue::INT_LIST)
      {
        column_name = column_name.substr(4);
        String value; //IntList value;
        Sql::extractValue<String>(&value, stmt, i); //IntList
        value = value.chop(1);
        value = value.substr(1);
        std::vector<String> value_list;
        IntList il = ListUtils::create<int>(value, ',');
        subordinate.setMetaValue(column_name, il);
        continue;
      }
      else if (column_type == DataValue::DOUBLE_LIST)
      {
        column_name = column_name.substr(4);
        String value; //DoubleList value;
        Sql::extractValue<String>(&value, stmt, i); //DoubleList
        value = value.chop(1);
        value = value.substr(1);          
        DoubleList dl = ListUtils::create<double>(value, ',');
        subordinate.setMetaValue(column_name, dl);
        continue;
      } 
      else if (column_type == DataValue::EMPTY_VALUE)
      {
        String value;
        Sql::extractValue<String>(&value, stmt, i);
        continue;
      }
    }


    // current_subordinate add bbox
    column_nr = cols_features + cols_subordinates + 1; //+ 1 = REF_ID
    ConvexHull2D hull = readBBox_(stmt, column_nr);
    subordinate.getConvexHulls().push_back(hull);

<<<<<<< HEAD
    std::cout << "subordinate id " << id_sub  << std::endl;
=======
    std::cout << "subordinate id " << id_sub << std::endl;
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
  
    return subordinate;
  }

/*
  // 
  std::pair<String, DataValue::DataType> getDataProcessing(sqlite3_stmt * stmt, int i)
  {
    String column_name = sqlite3_column_name(stmt, i);
    column_type = getColumnDatatype(column_name);

    if (column_type == DataValue::STRING_VALUE)
    {
      String value;
      column_name = column_name.substr(3);
      Sql::extractValue<String>(&value, stmt, i);
      std::make_pair(column_name, value); 
    } 
    else if (column_type == DataValue::INT_VALUE)
    {
      int value = 0;
      column_name = column_name.substr(3);
      Sql::extractValue<int>(&value, stmt, i);
      std::make_pair(column_name, value); 
    }
    else if (column_type == DataValue::DOUBLE_VALUE)
    {
      double value = 0.0;
      column_name = column_name.substr(3);
      Sql::extractValue<double>(&value, stmt, i);          
      std::make_pair(column_name, value); 
    }
    else if (column_type == DataValue::STRING_LIST)
    {
      String value; 
      Sql::extractValue<String>(&value, stmt, i);
      column_name = column_name.substr(4);
      StringList sl;
      // cut off "[" and "]""
      value = value.chop(1);
      value = value.substr(1);
      value.split(", ", sl);
      std::make_pair(column_name, value); 
    }
    else if (column_type == DataValue::INT_LIST)
    {
      String value; //IntList value;
      Sql::extractValue<String>(&value, stmt, i); //IntList
      column_name = column_name.substr(4);
      value = value.chop(1);
      value = value.substr(1);
      std::vector<String> value_list;
      IntList il = ListUtils::create<int>(value, ',');
      std::make_pair(column_name, value); 
    }
    else if (column_type == DataValue::DOUBLE_LIST)
    {
      String value; //DoubleList value;
      Sql::extractValue<String>(&value, stmt, i); //DoubleList
      column_name = column_name.substr(4);
      value = value.chop(1);
      value = value.substr(1);          
      DoubleList dl = ListUtils::create<double>(value, ',');
      std::make_pair(column_name, value); 
    }
    else if (column_type == DataValue::EMPTY_VALUE)
    {
      String value;
      column_name = column_name.substr(3);
      Sql::extractValue<String>(&value, stmt, i);
      std::make_pair(column_name, value); 
    }
  }
*/




// alternate version with table join 
/*
    // access of features and boundingboxes
    SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_TABLE_BOUNDINGBOX ON FEATURES_TABLE.id = FEATURES_TABLE_BOUNDINGBOX.ref_id
    // access of features, subordinates and boundingboxes
    SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_SUBORDINATES ON FEATURES_TABLE.id = FEATURES_SUBORDINATES.ref_id LEFT JOIN SUBORDINATES_TABLE_BOUNDINGBOX ON FEATURES_SUBORDINATES.id = SUBORDINATES_TABLE_BOUNDINGBOX.id

  // working modell
  read number of columns for all tables except dataprocessing
  combine access to joined tables via compound index
  one big loop with two different statements
  1. for features access
  2. for subordinate access
*/

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                   read FeatureMap as SQL database                                                    //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // read SQL database and store as FeatureMap
  // fitted snippet of FeatureXMLFile::load
  // based on TransitionPQPFile::readPQPInput

  FeatureMap FeatureSQLFile::read(const std::string& filename) const
  {
    //create empty FeatureMap object
    FeatureMap feature_map;

    sqlite3 *db;
    sqlite3_stmt * stmt;
    std::string select_sql;

    // Open database
    SqliteConnector conn(filename);
    db = conn.getDB();

    // set switches to access only existent tables
    bool features_switch = SqliteConnector::tableExists(db, "FEATURES_TABLE");
    bool subordinates_switch = SqliteConnector::tableExists(db, "FEATURES_SUBORDINATES");
    bool dataprocessing_switch = SqliteConnector::tableExists(db, "FEATURES_DATAPROCESSING");
    //bool features_bbox_switch = SqliteConnector::tableExists(db, "FEATURES_TABLE_BOUNDINGBOX");
    //bool subordinates_bbox_switch = SqliteConnector::tableExists(db, "SUBORDINATES_TABLE_BOUNDINGBOX");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                     store sqlite3 database as FeatureMap                                             //
    //                                                     1. features                                                                      //                                    
    //                                                     2. subordinates                                                                  //
    //                                                     3. boundingbox                                                                   //
    //                                                     4. dataprocessing                                                                //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    String sql;

    // get number of columns for each table except dataprocessing
    // 1. features
    sql = "SELECT * FROM FEATURES_TABLE;";
    SqliteConnector::prepareStatement(db, &stmt, sql);
    //SqliteConnector::executePreparedStatement(db, &stmt, sql);
    sqlite3_step(stmt);
    int cols_features = sqlite3_column_count(stmt);
    std::cout << "\n" << "Number of columns = " << cols_features << std::endl;
    sqlite3_finalize(stmt);

    // 2. subordinates
    sql = "SELECT * FROM FEATURES_SUBORDINATES;";
    SqliteConnector::prepareStatement(db, &stmt, sql);
    //SqliteConnector::executePreparedStatement(db, &stmt, sql);
    sqlite3_step(stmt);
    int cols_subordinates = sqlite3_column_count(stmt);
    std::cout << "Number of columns = " << cols_subordinates << std::endl;
    sqlite3_finalize(stmt);

    // 3. features bbox
    sql = "SELECT * FROM FEATURES_TABLE_BOUNDINGBOX;";
    SqliteConnector::prepareStatement(db, &stmt, sql);
    //SqliteConnector::executePreparedStatement(db, &stmt, sql);
    sqlite3_step(stmt);
    int cols_features_bbox = sqlite3_column_count(stmt);
    std::cout << "Number of columns = " << cols_features_bbox << std::endl;
    sqlite3_finalize(stmt);

    // 4. subordinate bbox
    sql = "SELECT * FROM SUBORDINATES_TABLE_BOUNDINGBOX;";
    SqliteConnector::prepareStatement(db, &stmt, sql);
    //SqliteConnector::executePreparedStatement(db, &stmt, sql);
    sqlite3_step(stmt);
    int cols_subordinates_bbox = sqlite3_column_count(stmt);
    std::cout << "Number of columns = " << cols_subordinates_bbox << std::endl;
    sqlite3_finalize(stmt);

    // 5. concatenation features bbox
    sql = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_TABLE_BOUNDINGBOX ON FEATURES_TABLE.id = FEATURES_TABLE_BOUNDINGBOX.ref_id;";
    SqliteConnector::prepareStatement(db, &stmt, sql);
    //SqliteConnector::executePreparedStatement(db, &stmt, sql);
    sqlite3_step(stmt);
    int cols_features_join_bbox = sqlite3_column_count(stmt);
    std::cout << "Number of columns = (features+bbox) " << cols_features_join_bbox << std::endl;
    sqlite3_finalize(stmt);

    // 6. concatenation features subordinates bbox 
    sql = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_SUBORDINATES ON FEATURES_TABLE.id = FEATURES_SUBORDINATES.ref_id LEFT JOIN SUBORDINATES_TABLE_BOUNDINGBOX ON FEATURES_SUBORDINATES.id = SUBORDINATES_TABLE_BOUNDINGBOX.id ORDER BY FEATURES_TABLE.ID;";
    SqliteConnector::prepareStatement(db, &stmt, sql);
    //SqliteConnector::executePreparedStatement(db, &stmt, sql);
    sqlite3_step(stmt);
    int cols_features_join_subordinates_join_bbox = sqlite3_column_count(stmt);
    std::cout << "Number of columns = (Features + Subordinates + bbox) " << cols_features_join_subordinates_join_bbox << std::endl;
    sqlite3_finalize(stmt);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                     store sqlite3 database as FeatureMap                                             //
    // traverse queries                                                                                                                     //
    // set all parameters: features, subordinates, boundingboxes                                                                            //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // add dataprocessing entries
    if (dataprocessing_switch)
    {
      /// 1. get feature data from database
      String sql = "SELECT * FROM  FEATURES_DATAPROCESSING;";
      SqliteConnector::prepareStatement(db, &stmt, sql);
      //SqliteConnector::executePreparedStatement(db, &stmt, sql);
      sqlite3_step(stmt);

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {

        // set dataprocessing parameters
        String id_s;
        size_t id = 0;

        // extract as String TODO
        Sql::extractValue<String>(&id_s, stmt, 0);
        std::istringstream iss(id_s);
        iss >> id;





        String software;
        Sql::extractValue<String>(&software, stmt, 1);
        String software_name;
        Sql::extractValue<String>(&software_name, stmt, 2);
        String data;
        Sql::extractValue<String>(&data, stmt, 3);
        String time;
        Sql::extractValue<String>(&time, stmt, 4);
        String action;
        Sql::extractValue<String>(&action, stmt, 5);

        // get values
        // id, RT, MZ, Intensity, Charge, Quality
        // save SQL column elements as feature
        feature_map.setUniqueId(id);
        DataProcessing dp;

        //software
        dp.getSoftware().setName(software);
        dp.getSoftware().setVersion(software_name);

        //time
        DateTime date_time;
        date_time.set(data + " " + time);
        dp.setCompletionTime(date_time);

        //actions
        std::set<DataProcessing::ProcessingAction> proc_actions;

        for (int i = 0; i < DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
        {
          if (action == DataProcessing::NamesOfProcessingAction[i])
          {
            proc_actions.insert((DataProcessing::ProcessingAction) i);
          }
        }
        dp.setProcessingActions(proc_actions);

        // access number of columns and infer DataType
        int cols = sqlite3_column_count(stmt);
        int column_type = 0;

        String column_name;

        for (int i = 5; i < cols ; ++i)
        {
          String column_name = sqlite3_column_name(stmt, i);
          column_type = getColumnDatatype(column_name);

          if (column_type == DataValue::STRING_VALUE)
          {
            String value;
            column_name = column_name.substr(3);
            Sql::extractValue<String>(&value, stmt, i);
            dp.setMetaValue(column_name, value); 
            continue;            
          } 
          else if (column_type == DataValue::INT_VALUE)
          {
            int value = 0;
            column_name = column_name.substr(3);
            Sql::extractValue<int>(&value, stmt, i);
            dp.setMetaValue(column_name, value); 
            continue;            
          }
          else if (column_type == DataValue::DOUBLE_VALUE)
          {
            double value = 0.0;
            column_name = column_name.substr(3);
            Sql::extractValue<double>(&value, stmt, i);          
            dp.setMetaValue(column_name, value); 
            continue;            
          }
          else if (column_type == DataValue::STRING_LIST)
          {
            String value; 
            Sql::extractValue<String>(&value, stmt, i);
            column_name = column_name.substr(4);
            StringList sl;
            // cut off "[" and "]""
            value = value.chop(1);
            value = value.substr(1);
            value.split(", ", sl);
            dp.setMetaValue(column_name, value); 
            continue;            
          }
          else if (column_type == DataValue::INT_LIST)
          {
            String value; //IntList value;
            Sql::extractValue<String>(&value, stmt, i); //IntList
            column_name = column_name.substr(4);
            value = value.chop(1);
            value = value.substr(1);
            std::vector<String> value_list;
            IntList il = ListUtils::create<int>(value, ',');
            dp.setMetaValue(column_name, value); 
            continue;            
          }
          else if (column_type == DataValue::DOUBLE_LIST)
          {
            String value; //DoubleList value;
            Sql::extractValue<String>(&value, stmt, i); //DoubleList
            column_name = column_name.substr(4);
            value = value.chop(1);
            value = value.substr(1);          
            DoubleList dl = ListUtils::create<double>(value, ',');
            dp.setMetaValue(column_name, value); 
            continue;            
          }
          else if (column_type == DataValue::EMPTY_VALUE)
          {
            String value;
            column_name = column_name.substr(3);
            Sql::extractValue<String>(&value, stmt, i);
            dp.setMetaValue(column_name, value); 
            continue;            
          }
        }

        feature_map.getDataProcessing().push_back(dp);
        sqlite3_step(stmt);

      }

      sqlite3_finalize(stmt);

    } //dataprocessing_switch








    std::map<int64_t, size_t> map_fid_to_index; 
  
    // add feature entries
    if (features_switch)
    {
      /// 1. get feature data from database
      String sql = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_TABLE_BOUNDINGBOX ON FEATURES_TABLE.id = FEATURES_TABLE_BOUNDINGBOX.ref_id;";
      SqliteConnector::prepareStatement(db, &stmt, sql);
      //SqliteConnector::executePreparedStatement(db, &stmt, sql);
      sqlite3_step(stmt);
      
      Feature* feature = nullptr;

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        int bbox_col = cols_features + 1 + 4; //56 in features + 1 offset id + 4 position of bbox_idx column
        int bbox_idx = 0;
        Sql::extractValue<int>(&bbox_idx, stmt, bbox_col);

        if (bbox_idx == 0)
        {
          feature_map.push_back(Feature());
          feature = &feature_map[feature_map.size() - 1];

          // set feature parameters
          String id_s;
          int64_t id = 0;

          // extract as String TODO
          Sql::extractValue<String>(&id_s, stmt, 0);
          std::istringstream iss(id_s);
          iss >> id;

          double rt = 0.0;
          Sql::extractValue<double>(&rt, stmt, 1);
          double mz = 0.0;
          Sql::extractValue<double>(&mz, stmt, 2);
          double intensity = 0.0;
          Sql::extractValue<double>(&intensity, stmt, 3);
          int charge = 0;
          Sql::extractValue<int>(&charge, stmt, 4);
          double quality = 0.0;
          Sql::extractValue<double>(&quality, stmt, 5);

          // get values id, RT, MZ, Intensity, Charge, Quality
          // save SQL column elements as feature
          feature->setUniqueId(id);
          feature->setRT(rt);
          feature->setMZ(mz);
          feature->setIntensity(intensity);
          feature->setCharge(charge);
          feature->setOverallQuality(quality);

          // index to account for offsets in joined tables
          int col_nr = 6;
          // save userparam columns
          for (int i = col_nr; i < cols_features ; ++i)
          {
            String column_name = sqlite3_column_name(stmt, i);
            int column_type = getColumnDatatype(column_name);

            //std::cout << "\n" << column_name << std::endl;

            if (column_type == DataValue::STRING_VALUE)
            {
              column_name = column_name.substr(3);
              String value;
              Sql::extractValue<String>(&value, stmt, i);
              feature->setMetaValue(column_name, value); 
              continue;
            } 
            else if (column_type == DataValue::INT_VALUE)
            {
              column_name = column_name.substr(3);
              int value = 0;
              Sql::extractValue<int>(&value, stmt, i);
              feature->setMetaValue(column_name, value); 
              continue;
            } 
            else if (column_type == DataValue::DOUBLE_VALUE)
            {
              column_name = column_name.substr(3);
              double value = 0.0;
              Sql::extractValue<double>(&value, stmt, i);          
              feature->setMetaValue(column_name, value); 
              continue;
            } 
            else if (column_type == DataValue::STRING_LIST)
            {
              column_name = column_name.substr(4);
              String value; 
              Sql::extractValue<String>(&value, stmt, i);

              StringList sl;
              // cut off "[" and "]""
              value = value.chop(1);
              value = value.substr(1);
              value.split(", ", sl);
              feature->setMetaValue(column_name, sl);
              continue;
            } 
            else if (column_type == DataValue::INT_LIST)
            {
              column_name = column_name.substr(4);
              String value; //IntList value;
              Sql::extractValue<String>(&value, stmt, i); //IntList
              value = value.chop(1);
              value = value.substr(1);
              std::vector<String> value_list;
              IntList il = ListUtils::create<int>(value, ',');
              feature->setMetaValue(column_name, il);
              continue;
            } 
            else if (column_type == DataValue::DOUBLE_LIST)
            {
              column_name = column_name.substr(4);
              String value; //DoubleList value;
              Sql::extractValue<String>(&value, stmt, i); //DoubleList
              value = value.chop(1);
              value = value.substr(1);          
              DoubleList dl = ListUtils::create<double>(value, ',');
              feature->setMetaValue(column_name, dl);
              continue;
            } 
            else if (column_type == DataValue::EMPTY_VALUE)
            {
              String value;
              Sql::extractValue<String>(&value, stmt, i);
              continue;
            }
          }
        
          int column_nr = cols_features;
          ConvexHull2D hull = readBBox_(stmt, column_nr);
          feature->getConvexHulls().push_back(hull);

          // map feature id to index in feature map (needed to map subordinates to features)
          cout << "Adding feature ID with map index: " << feature->getUniqueId() << "\t" << feature_map.size() - 1 << endl;
          map_fid_to_index[feature->getUniqueId()] = feature_map.size() - 1; 
        } 
        else if (bbox_idx > 0) // add new bb to existing feature
        {
          std::cout << "\n bbox_idx = " << bbox_idx << std::endl;

          int column_nr = cols_features;
          ConvexHull2D hull = readBBox_(stmt, column_nr);
          feature->getConvexHulls().push_back(hull);
        }
        sqlite3_step(stmt);
<<<<<<< HEAD

=======
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
      }
      sqlite3_finalize(stmt);
    }

<<<<<<< HEAD
    // add subordinate entries
=======
    // subordinates
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
    if (subordinates_switch)
    {
      // join features with subordinates + subordinate bounding boxes
      sql = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_SUBORDINATES ON FEATURES_TABLE.id = FEATURES_SUBORDINATES.ref_id LEFT JOIN SUBORDINATES_TABLE_BOUNDINGBOX ON FEATURES_SUBORDINATES.id = SUBORDINATES_TABLE_BOUNDINGBOX.id ORDER BY FEATURES_TABLE.ID;";
      SqliteConnector::prepareStatement(db, &stmt, sql);
      sqlite3_step(stmt);

      int column_nr = 0;
<<<<<<< HEAD
      int64_t f_id_prev = 0;
=======
      size_t f_id_prev = 0;
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        // get feature ID
        int64_t f_id = 0;
        String f_id_string;
        Sql::extractValue<String>(&f_id_string, stmt, 0);
        std::istringstream iss(f_id_string);
        iss >> f_id;

        // get boundingbox index
        int bbox_col = cols_features_join_subordinates_join_bbox - 1; // - 1 to address last column
        int bbox_idx = 0; // set start idx of boundingbox
        Sql::extractValue<int>(&bbox_idx, stmt, bbox_col);


        Feature* feature = &feature_map[map_fid_to_index[f_id]];
<<<<<<< HEAD

        std::vector<Feature>* subordinates = &feature->getSubordinates();

=======
        std::vector<Feature>* subordinates = &feature->getSubordinates();
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5

        if (f_id != f_id_prev) // new feature
        {
          f_id_prev = f_id;
          column_nr = cols_features + 1 + 1;// size features + column SUB_IDX + column REF_ID
          subordinates->push_back(readSubordinate_(stmt, column_nr, cols_features, cols_subordinates)); // add subordinate and first bounding box (if present)
<<<<<<< HEAD
          sqlite3_step(stmt);
          continue;
=======
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
        }
        else // same feature, different subordinate or bounding box
        {

          if (bbox_idx == 0) 
          {
            column_nr = cols_features + 1 + 1;// size features + column SUB_IDX + column REF_ID
            subordinates->push_back(readSubordinate_(stmt, column_nr, cols_features, cols_subordinates)); // read subordinate + first bounding box
          }
          else // add new bounding box to current subordinate
          {
            column_nr = cols_features + cols_subordinates + 1; //+ 1 = REF_ID
            (*subordinates)[subordinates->size()-1].getConvexHulls().push_back(readBBox_(stmt, column_nr));
          }
<<<<<<< HEAD
          sqlite3_step(stmt);       
          continue;
=======
>>>>>>> 1ea84ab1ac4adad0108185a099238e3fe621a0c5
        }
        sqlite3_step(stmt);
      } // while
      sqlite3_finalize(stmt);
    }  // subordinate switch

    return feature_map;
  }  // end of FeatureSQLFile::read
} // namespace OpenMS




































//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                     store sqlite3 database as FeatureMap                                             //
    /// alterantive version where 2 statements are used to access the joined tables
    // and traversed  once
    // read features and bounding boxes
    // read subordinates and bounding boxes
    // link features + bbox and subordiates + bbox
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // pseudo code box 
    /*
      feature
        iterate over identical features and subordinates
        {
          if f_id == f_id_prev
          {
            if s_id == s_id_prev 
            {
              if s_bbox_idx == s_bbox_idx_prev
              {
                push back bbox to current subordinate
              }
              else 
              {
                f_id same, s_id same, bbox same -> must be same entry -> is impossible bc table conflict
                use as catch case for NULL s_bbox entries 
              }
            }
            else if (s_id != s_id_prev)
            {
              push back subordinate to subordinates
              clear bbox vector since new subordinate = new convexhulls
              
              if ((s_bbox_idx == s_bbox_idx_prev) && (s_bbox_idx != 0))
              {
                use as catch case for NULL s_bbox entries
              }
              else ((s_bbox_idx != s_bbox_idx_prev) || (s_bbox_idx == 0))
              {
                push_back s_bbox to bbox vector
              }            
            }  
          }
          else if (f_id != f_id_prev)
          {
            push back feature entries to feature_map
            feature = Feature();

            get bbox since test for subordinates is superfluous, new f_id -> new subordinate
            0. clear feature, clear subordinate, clear subordinates vector, clear bbox vector
            // feature = Feature();
            // subordinate = Feature();
            // subordinates.clear();
            // bbox_vector.clear(); ? subordinates and bbox one object?
            1. get f_bbox
            2. get subordinate + f_bbox 
            3. set f_id_prev = f_id, s_id_prev = s_id 

            save to featuremap

          }

          step f_line
          step s_line
        }
    */
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/*

    if (features_switch)
    {
      sqlite3 *db;
      sqlite3_stmt * stmt_features;
      sqlite3_stmt * stmt_subordinates;
      
      // Open database
      SqliteConnector conn(filename);
      db = conn.getDB();

      /// 1. get feature data from database
      // unordered
      //String sql_feats = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_TABLE_BOUNDINGBOX ON FEATURES_TABLE.id = FEATURES_TABLE_BOUNDINGBOX.ref_id;";
      String sql_feats = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_TABLE_BOUNDINGBOX ON FEATURES_TABLE.id = FEATURES_TABLE_BOUNDINGBOX.ref_id ORDER BY FEATURES_TABLE.id;";
      String sql_subs = "SELECT * FROM  FEATURES_TABLE LEFT JOIN FEATURES_SUBORDINATES ON FEATURES_TABLE.id = FEATURES_SUBORDINATES.ref_id LEFT JOIN SUBORDINATES_TABLE_BOUNDINGBOX ON FEATURES_SUBORDINATES.id = SUBORDINATES_TABLE_BOUNDINGBOX.id ORDER BY FEATURES_TABLE.ID;";

      SqliteConnector::prepareStatement(db, &stmt_features, sql_feats);
      SqliteConnector::prepareStatement(db, &stmt_subordinates, sql_subs);
     
      sqlite3_step(stmt_features);
      sqlite3_step(stmt_subordinates);

      // initialization block
      // feature and subordinate stmts used to push_back columns of both tables
      // set starting values
      Feature feature;
      Feature subordinate;
      std::vector<Feature> subordinates;
      // ids
      long f_id_prev = 0;
      long s_id_prev = 0;
      // subordinate index
      int sub_idx_prev = 0;
      // bboxes
      int f_bbox_idx_prev = 0;
      int s_bbox_idx_prev = 0;

      // first iteration variable to catch enter loop case
      // use loop_switch
      // set to false (= not entered) on start
      // and to false if last line is reached
      bool first_feat_switch = false;

      // TODO set last_feature_switch ? ###########################

      // provisional code
      // test counter to access only first 2 lines in joined tables
      int test_counter = 1;

      // get features + bboxes, subordinates + bboxes in a single loop
      for (int i = 0; i <= test_counter; ++i)
      {
        std::cout << "\n" << i << std::endl;

      

        /// get data starts here
        // ##################################################
        // get data for features and bbox
        // get feature id
        String feat_id;
        long f_id = 0;
        Sql::extractValue<String>(&feat_id, stmt_features, 0);
        f_id = stol(feat_id);

        // get feature bbox_idx
        int f_bbox_col = cols_features + 1 + 4; //56 in features + 1 offset id + 4 position of bbox_idx column
        int f_bbox_idx = 0;
        Sql::extractValue<int>(&f_bbox_idx, stmt_features, f_bbox_col);

        // ##################################################
        // get data for subordinates and bbox
        // get subordinate ID
        String sub_id;
        long s_id = 0;
        Sql::extractValue<String>(&sub_id, stmt_subordinates, cols_features);
        s_id = stol(sub_id);        

        // return sub_idx (n-th subordinate)
        int sub_idx = 0;
        Sql::extractValue<int>(&sub_idx, stmt_subordinates, 57);

        // get subordinate bbox_idx (n-th bbox)
        int s_bbox_col = cols_features_join_subordinates_join_bbox - 1; // - 1 to address last column
        int s_bbox_idx = 0; // set start idx of boundingbox
        Sql::extractValue<int>(&s_bbox_idx, stmt_subordinates, s_bbox_col);

        // ##################################################
        // first run case
        if (first_feat_switch == false)
          {
            first_feat_switch = true;
            f_id_prev = f_id;
            s_id_prev = s_id;

            sub_idx_prev = sub_idx;

            f_bbox_idx_prev = f_bbox_idx;
            s_bbox_idx_prev = s_bbox_idx;
          }

        // ##################################################
        // iterate over identical features and subordinates
        // ##################################################

        std::cout << "f_id, f_id_prev, s_id, s_id_prev " << f_id << "\t" << f_id_prev << "\t" << s_id << "\t" << s_id_prev << std::endl;

        if (f_id == f_id_prev) 
        {
          std::cout << "f_id == f_id_prev" << std::endl;


          if (s_id == s_id_prev)
          {
            std::cout << "s_id == s_id_prev" << std::endl;
            if (s_bbox_idx == s_bbox_idx_prev)
            {
              //push back bbox to current subordinate
              std::cout << "s_bbox_idx == s_bbox_idx_prev" << std::endl;
            }
            else if (s_bbox_idx != s_bbox_idx_prev) 
            {
              std::cout << "f_id same, s_id same, bbox same -> must be same entry -> is impossible bc table conflict" << std::endl;
            }
          }
          else if (s_id != s_id_prev)
          {
            std::cout << "s_id != s_id_prev" << std::endl;
            // push back subordinate to subordinates
            subordinates.push_back(subordinate);
            // clear bbox vector since new subordinate = new convexhulls
            subordinates.clear();
            
            if ((s_bbox_idx == s_bbox_idx_prev) && (s_bbox_idx != 0))
            {
              // use as catch case for NULL s_bbox entries
              std::cout << "if ((s_bbox_idx == s_bbox_idx_prev) && (s_bbox_idx != 0)) " << std::endl;
            }
            else if ((s_bbox_idx != s_bbox_idx_prev) || (s_bbox_idx == 0))
            {
              // push_back s_bbox to bbox vector
              std::cout << "subordinate.getConvexHulls().push_back(hull);" << std::endl;
            }            
          }  
        }
        else if (f_id != f_id_prev)
        {
          std::cout << "f_id != f_id_prev" << std::endl;
        
          // push back feature entries to feature_map
          // feature = Feature();

          // get bbox since test for subordinates is superfluous, new f_id -> new subordinate
          // 0. clear feature, clear subordinate, clear subordinates vector, clear bbox vector
          // // feature = Feature();
          // // subordinate = Feature();
          // // subordinates.clear();
          // // bbox_vector.clear(); ? subordinates and bbox one object?
          // 1. get f_bbox
          // 2. get subordinate + f_bbox 
          // 3. set f_id_prev = f_id, s_id_prev = s_id 

          // save to featuremap
        
        }

        f_id_prev = f_id;
        s_id_prev = s_id;

        sub_idx_prev = sub_idx;

        f_bbox_idx_prev = f_bbox_idx;
        s_bbox_idx_prev = s_bbox_idx;

        //step f_line
        sqlite3_step(stmt_features);
        //step s_line
        sqlite3_step(stmt_subordinates);

      } // close loop

      // closing statements for features and subordinates query
      sqlite3_finalize(stmt_features);
      sqlite3_finalize(stmt_subordinates);

    } // close features_switch

    return feature_map;

  } // end of FeatureSQLFile::read

} // namespace OpenMS


*/

