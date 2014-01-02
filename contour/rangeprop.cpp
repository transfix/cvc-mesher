
// $Id: rangeprop.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <stdio.h>
#include <string.h>
#include "rangeprop.h"

#define DEBUGNo

extern int verbose;

void
rangeProp::compSeeds(void)
{
   RangePropRec rpr, *item, *qitem;
   Range fullrange, added, outgoing;
   int index;
   float min, max;
   int adjc;
   int done;
   u_int c;

   if (verbose)
   printf("------- computing seeds\n");

   // clear the array of mark bits
   plot.ClearTouched();
   seeds.Clear();

   // insert cell 0 into queue to begin
   rpr.cellid = 0;
   data.getCellRange(0, min, max);
   rpr.resp.Set(min, max);
   rpr.comp.MakeEmpty();
   queue.enqueue(rpr);

   done = 0;

   // process queue of cells
   while (! queue.isEmpty() ) {
      // get the item
      item = queue.dequeue();

      done++;

      if (verbose)
      if (done % 1000 == 0)
         printf("%g%% done\n", 100.0*done/(float)data.getNCells());
#ifdef DEBUG
{ int i;
sleep(1);
printf("processing cell %d\n", item->cellid);
printf(" verts: ");
for (i=0; i<data.getNCellVerts(); i++)
   printf("%d ", data.getCellVert(item->cellid, i));
printf(" values: ");
for (i=0; i<data.getNCellVerts(); i++)
   printf("%f ", data.getValue(data.getCellVert(item->cellid, i)));
printf("\nresp: ");
item->resp.Print();
printf("comp: ");
item->comp.Print(); }
#endif


      // mark this cell as processed
      plot.TouchCell(item->cellid);

#ifdef DEBUG
printf("touched cell %d\n", item->cellid);
#endif

      // compute the outgoing range (range which can be further propagated)
      outgoing.MakeEmpty();
#ifdef DEBUG
printf("%d adjacencies\n", data.getNCellFaces());
#endif


      for (c=0; c<data.getNCellFaces(); c++) {
         adjc = data.getCellAdj(item->cellid, c);

#ifdef DEBUG
printf("adj cell %d\n", adjc);
#endif

         if (adjc != -1 && !plot.CellTouched(adjc)) {
            // the range of the shared face may be propagated
            data.getFaceRange(item->cellid, c, min, max);
            outgoing += Range(min,max);
         }
      }

#ifdef DEBUG
printf("outgoing range: ");
outgoing.Print();
#endif

      // this is the full range of responsibility
      fullrange = (item->resp + outgoing) - item->comp;

#ifdef DEBUG
printf("fullrange: ");
fullrange.Print();
#endif

      if (fullrange.Empty() || (!outgoing.Empty() &&
                             (outgoing.MinAll() <= fullrange.MinAll() &&
                              outgoing.MaxAll() >= fullrange.MaxAll()))) {
#ifdef DEBUG
printf("propagated\n");
#endif
         // propagate entire range
         for (c=0; c<data.getNCellFaces(); c++) {
            adjc = data.getCellAdj(item->cellid, c);

            if (adjc != -1 && !plot.CellTouched(adjc)) {
               data.getFaceRange(item->cellid, c, min, max);

               // compute the range which should be propagated to this cell
               added = Range(min,max) - item->comp;
               item->comp += Range(min,max);

#ifdef DEBUG
printf("   adj %d:\n", adjc);
#endif

               rpr.cellid = adjc;
               if ((index=queue.find(rpr)) != -1) {
                  qitem=queue.getItem(index);
                  qitem->resp += added;
                  qitem->comp = Range(min,max) - added;
#ifdef DEBUG
printf("      found cell, new resp: ");
qitem->resp.Print();
printf("      new comp: ");
qitem->comp.Print();
#endif
               }
               else {
                  rpr.resp = added;
                  rpr.comp = Range(min,max) - added;
                  queue.enqueue(rpr);
#ifdef DEBUG
printf("      new cell, resp: ");
rpr.resp.Print();
printf("      comp: ");
rpr.comp.Print();
#endif
               }
            }
         }
      }
      else {
#ifdef DEBUG
printf("seed cell\n");
#endif
         // cell is a seed cell
         seeds.AddSeed(item->cellid, fullrange.MinAll(), fullrange.MaxAll());

         // set range of faces to the _complement_ of adjacent cells
         for (c=0; c<data.getNCellFaces(); c++) {
            adjc = data.getCellAdj(item->cellid, c);

            if (adjc != -1 && !plot.CellTouched(adjc)) {
               data.getFaceRange(item->cellid, c, min, max);
               rpr.cellid = adjc;
               if ((index=queue.find(rpr)) != -1) {
                  qitem=queue.getItem(index);
                  qitem->comp.Set(min,max);
//                  qitem->comp += Range(min,max);
// note : can we ADD to the comp here?
               }
               else {
                  rpr.resp.MakeEmpty();
                  rpr.comp.Set(min,max);
                  queue.enqueue(rpr);
               }
            }
         }
      }
   }

//   for (c=0; c<data.getNCells(); c++)
//      if (!plot.CellTouched(c))
//         printf("cell %d not touched\n", c);

   if (verbose)
   printf("computed %d seeds\n", seeds.getNCells());
}
