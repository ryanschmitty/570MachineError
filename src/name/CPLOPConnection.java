package name;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.sql.SQLException;

/**
 * A class for accessing the CPLOP Databse.
 * The connection is persisted over calls and is
 *  automatically reestablished if it becomes stale.
 */
public class CPLOPConnection
{
   private static final String DB_DRIVER = "com.mysql.jdbc.Driver";
   private static final String DB_URL = "jdbc:mysql://abra.csc.calpoly.edu/cplop?autoReconnect=true";
   private static final String DB_USER = "csc570";
   private static final String DB_PASS = "ilovedata570";

   private Connection conn;

   /**
    * Test the connection.
    */
   public static void main(String[] args)
   {
      try
      {
         CPLOPConnection cplop = new CPLOPConnection();

         List<Map<String, Object>> res;
         
         /*
         res = cplop.getIsolates();
         for (Map<String, Object> isolate : res)
         {
            System.out.println("Isolate: " + isolate.get("isolate") + 
             ", count: " + isolate.get("count"));
         }
         */
         
         /*
         res = cplop.getPyroprints("Hu-1526");
         for (Map<String, Object> pyroprint : res)
         {
            System.out.println(pyroprint);
         }
         */
         
         res = cplop.getHistogram(2377);
         for (Map<String, Object> pyroprint : res)
         {
            System.out.println(pyroprint);
         }
      }
      catch (Exception ex)
      {
         ex.printStackTrace();
      }
   }

   /**
    * @throws SQLException if the DriverManager can't get a connection.
    * @throws DriverException if there is a problem instantiating the DB driver.
    */
   public CPLOPConnection() throws SQLException, DriverException
   {
      try
      {
         Class.forName(DB_DRIVER);

         conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);
      }
      catch (ClassNotFoundException classEx)
      {
         throw new DriverException("Unable to instantiate DB Driver: " + DB_DRIVER);
      }
   }

   /**
    * Get all the different isolates and the number of 
    *  pyroprints that are in the database fop each.
    * Note that the number is not the number of times the isolate
    *  was pyroprinted, but the number of pyroprints that are in the DB.
    *
    * @return A list of hashes representing the isolates and their ids.
    *  Each row is a Map that maps attribute names to values.
    *  'isolate' maps to the isolate id (String).
    *  'count' maps to the number of pyroprints for the isolate (Integer).
    *
    * @throws SQLException if the query fails.
    */
   public List<Map<String, Object>> getIsolates() throws SQLException
   {
      List<Map<String, Object>> rtn = new ArrayList<Map<String, Object>>();
      Statement stmt = null;
      ResultSet results = null;

      String query = "SELECT isoID, COUNT(*)" + 
       " FROM Isolates i JOIN Pyroprints p USING(isoID)" + 
       " WHERE isoID IN (SELECT DISTINCT(isoID) FROM Histograms)" +
       " GROUP BY isoID" +
       " ORDER BY count(*) DESC, isoID";

      try
      {
         stmt = conn.createStatement();

         results = stmt.executeQuery(query);

         while (results.next())
         {
            Map<String, Object> tuple = new HashMap<String, Object>();
            tuple.put("isolate", results.getString(1)); 
            tuple.put("count", new Integer(results.getInt(2))); 

            rtn.add(tuple);
         }
      }
      catch (SQLException sqlEx)
      {
         //Rethrow the exception
         throw sqlEx;
      }
      finally
      {
         if (results != null)
         {
            results.close();
         }

         if (stmt != null)
         {
            stmt.close();
         }
      }

      return rtn;
   }


   /**
    * Given an isolate (isolate id), give the data for all pyroprints in the
    *  the databse.
    *
    * @param isoId The islate id. The same one given from getIsolate.
    *
    * @return A list of hashes representing the pyroprints.
    *  Each row is a Map that maps attribute names to values.
    *  'pyroprint' -> pyroprint's id (Integer).
    *  'region' -> unique identifier for the region. (String).
    *  'forwardPrimer' -> unique identifier for the forward primer (String).
    *  'reversePrimer' -> unique identifier for the reverse primer (String).
    *  'sequencePrimer' -> unique identifier for the sequence primer (String).
    *  'dispensation' -> unique identifier for the dispensation sequence (String).
    *
    * @throws SQLException if the query fails.
    */
   public List<Map<String, Object>> getPyroprints(String isoId) throws SQLException
   {
      List<Map<String, Object>> rtn = new ArrayList<Map<String, Object>>();
      Statement stmt = null;
      ResultSet results = null;

      String query = String.format(
       "SELECT pyroID, appliedRegion, dsName, forPrimer, revPrimer, seqPrimer" +
       " FROM Pyroprints" +
       " WHERE isoID = '%s'",
       isoId);

      try
      {
         stmt = conn.createStatement();

         results = stmt.executeQuery(query);

         while (results.next())
         {
            Map<String, Object> tuple = new HashMap<String, Object>();
            tuple.put("pyroprint", new Integer(results.getInt(1))); 
            tuple.put("region", results.getString(2)); 
            tuple.put("dispensation", results.getString(3)); 
            tuple.put("forwardPrimer", results.getString(4)); 
            tuple.put("reversePrimer", results.getString(5)); 
            tuple.put("sequencePrimer", results.getString(6)); 

            rtn.add(tuple);
         }
      }
      catch (SQLException sqlEx)
      {
         //Rethrow the exception
         throw sqlEx;
      }
      finally
      {
         if (results != null)
         {
            results.close();
         }

         if (stmt != null)
         {
            stmt.close();
         }
      }

      return rtn;
   }

   /**
    * Given a pyroprint (pyroprint id), give the histogram.
    *
    * @param pyroId The pyroprint id. The same one given from getPyroprints.
    *
    * @return A list of hashes representing the histogram.
    *  Each row is a Map that maps attribute names to values.
    *  'peakHeight' -> the peak height for this position. (Double)
    *  'peakArea' -> the peak area for this position. (Double)
    *  'peakWidth' -> the width of the peak for this position. (Double)
    *  'compensatedPeakHeight' -> the compensated peak height for this position. (Double)
    *  'nucleotide' - The nucleotide at this position. (String)
    *
    * @throws SQLException if the query fails.
    */
   public List<Map<String, Object>> getHistogram(int pyroId) throws SQLException
   {
      List<Map<String, Object>> rtn = new ArrayList<Map<String, Object>>();
      Statement stmt = null;
      ResultSet results = null;

      String query = String.format(
       "SELECT pHeight, PeakArea, PeakWidth, cPeakHeight, nucleotide" +
       " FROM Histograms" + 
       " WHERE pyroID = %d" +
       " ORDER BY position",
       pyroId);

      try
      {
         stmt = conn.createStatement();

         results = stmt.executeQuery(query);

         while (results.next())
         {
            Map<String, Object> tuple = new HashMap<String, Object>();
            tuple.put("peakHeight", new Double(results.getDouble(1)));
            tuple.put("peakArea", new Double(results.getDouble(2)));
            tuple.put("peakWidth", new Double(results.getDouble(3)));
            tuple.put("compensatedPeakHeight", new Double(results.getDouble(4)));
            tuple.put("nucleotide", results.getString(5));

            rtn.add(tuple);
         }
      }
      catch (SQLException sqlEx)
      {
         //Rethrow the exception
         throw sqlEx;
      }
      finally
      {
         if (results != null)
         {
            results.close();
         }

         if (stmt != null)
         {
            stmt.close();
         }
      }

      return rtn;
   }


   /**
    * Perform just a general query.
    * Make sure to close your result set after you are finished with it.
    *
    * @throws SQLException if the query generates an error.
    *
    * @return the ResultSet that came from the query.
    */
   public ResultSet generalQuery(String query) throws SQLException
   {
      Statement stmt = conn.createStatement();

      return stmt.executeQuery(query);
   }

   /**
    * Exception for when there is a problem with the DB driver.
    */
   public class DriverException extends Exception
   {
      public DriverException(String msg)
      {
         super(msg);
      }
   }
}
