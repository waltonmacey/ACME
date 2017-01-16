/**
 * @file
 * Algorithms modeled after spmd_utils in the Community
 * Atmosphere Model; C translation. This includes MPI_Gather,
 * MPI_Gatherv, and MPI_Alltoallw with flow control options.
 *
 * @author Jim Edwards
 * @date 2014
 */
#include <config.h>
#include <pio.h>
#include <pio_internal.h>

/**
 * Returns the smallest power of 2 greater than i.
 *
 * @param i input number
 * @returns the smallest power of 2 greater than i.
 */
int ceil2(int i)
{
    int p = 1;
    while (p < i)
        p *= 2;

    return(p);
}

/**
 * Given integers p and k between 0 and np-1,
 * if (p+1)^k <= np-1 then return (p+1)^k else -1.
 *
 * @param np
 * @param p integer between 0 and np - 1.
 * @param k integer between 0 and np - 1.
 * @returns (p + 1) ^ k else -1.
 */
int pair(int np, int p, int k)
{
    int q = (p + 1) ^ k ;
    int pair = (q > np - 1) ? -1 : q;
    return pair;
}

/**
 * Provides the functionality of MPI_Alltoallw with flow control
 * options. Generalized all-to-all communication allowing different
 * datatypes, counts, and displacements for each partner
 *
 * @param sendbuf starting address of send buffer
 * @param sendcounts integer array equal to the number of tasks in
 * communicator comm (ntasks). It specifies the number of elements to
 * send to each processor
 * @param sdispls integer array (of length ntasks). Entry j
 * specifies the displacement in bytes (relative to sendbuf) from
 * which to take the outgoing data destined for process j.
 * @param sendtypes array of datatypes (of length ntasks). Entry j
 * specifies the type of data to send to process j.
 * @param recvbuf address of receive buffer.
 * @param recvcounts integer array (of length ntasks) specifying the
 * number of elements that can be received from each processor.
 * @param rdispls integer array (of length ntasks). Entry i
 * specifies the displacement in bytes (relative to recvbuf) at which
 * to place the incoming data from process i.
 * @param recvtypes array of datatypes (of length ntasks). Entry i
 * specifies the type of data received from process i.
 * @param comm MPI communicator for the MPI_Alltoallw call.
 * @param handshake if true, use handshaking.
 * @param isend the isend bool indicates whether sends should be
 * posted using mpi_irsend which can be faster than blocking
 * sends. When flow control is used max_requests > 0 and the number of
 * irecvs posted from a given task will not exceed this value. On some
 * networks too many outstanding irecvs will cause a communications
 * bottleneck.
 * @param max_requests If 0, no flow control is used.
 * @returns 0 for success, error code otherwise.
 */
int pio_swapm(void *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype *sendtypes,
              void *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype *recvtypes,
              MPI_Comm comm, bool handshake, bool isend, int max_requests)
{
    int ntasks;  /* Number of tasks in communicator comm. */
    int my_rank; /* Rank of this task in comm. */
    int tag;
    int offset_t;
    int steps;
    int istep;
    int rstep;
    int p;
    int maxreq;
    int maxreqh;
    int hs = 1; /* Used for handshaking. */
    void *ptr;
    MPI_Status status; /* Not actually used - replace with MPI_STATUSES_IGNORE. */
    int mpierr;  /* Return code from MPI functions. */

    LOG((2, "pio_swapm handshake = %d isend = %d max_requests = %d", handshake,
         isend, max_requests));

    /* Get my rank and size of communicator. */
    if ((mpierr = MPI_Comm_size(comm, &ntasks)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    if ((mpierr = MPI_Comm_rank(comm, &my_rank)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);

    LOG((2, "ntasks = %d my_rank = %d", ntasks, my_rank));

    /* Now we know the size of these arrays. */
    int swapids[ntasks];
    MPI_Request rcvids[ntasks];
    MPI_Request sndids[ntasks];
    MPI_Request hs_rcvids[ntasks];

    /* Print some debugging info, if logging is enabled. */
#if PIO_ENABLE_LOGGING
    {
        for (int p = 0; p < ntasks; p++)
            LOG((3, "sendcounts[%d] = %d", p, sendcounts[p]));
        for (int p = 0; p < ntasks; p++)
            LOG((3, "sdispls[%d] = %d", p, sdispls[p]));
        for (int p = 0; p < ntasks; p++)
            LOG((3, "sendtypes[%d] = %d", p, sendtypes[p]));
        for (int p = 0; p < ntasks; p++)
            LOG((3, "recvcounts[%d] = %d", p, recvcounts[p]));
        for (int p = 0; p < ntasks; p++)
            LOG((3, "rdispls[%d] = %d", p, rdispls[p]));
        for (int p = 0; p < ntasks; p++)
            LOG((3, "recvtypes[%d] = %d", p, recvtypes[p]));
    }
#endif /* PIO_ENABLE_LOGGING */

    /* If max_requests == 0 no throttling is requested and the default
     * mpi_alltoallw function is used. */
    if (max_requests == 0)
    {
        /* Call the MPI alltoall without flow control. */
        LOG((3, "Calling MPI_Alltoallw without flow control."));
        if ((mpierr = MPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf,
                                    recvcounts, rdispls, recvtypes, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        return PIO_NOERR;
    }

    /* an index for communications tags */
    offset_t = ntasks;

    /* Send to self. */
    if (sendcounts[my_rank] > 0)
    {
        void *sptr, *rptr;
        tag = my_rank + offset_t;
        sptr = (char *)sendbuf + sdispls[my_rank];
        rptr = (char *)recvbuf + rdispls[my_rank];

        /*
          MPI_Type_get_extent(sendtypes[my_rank], &lb, &extent);
          printf("%s %d %d %d\n",__FILE__,__LINE__,extent, lb);
          MPI_Type_get_extent(recvtypes[my_rank], &lb, &extent);
          printf("%s %d %d %d\n",__FILE__,__LINE__,extent, lb);
        */

#ifdef ONEWAY
        /* If ONEWAY is true we will post mpi_sendrecv comms instead
         * of irecv/send. */
        if ((mpierr = MPI_Sendrecv(sptr, sendcounts[my_rank],sendtypes[my_rank],
                                   my_rank, tag, rptr, recvcounts[my_rank], recvtypes[my_rank],
                                   my_rank, tag, comm, &status)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
#else
        if ((mpierr = MPI_Irecv(rptr, recvcounts[my_rank], recvtypes[my_rank],
                                my_rank, tag, comm, rcvids)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Send(sptr, sendcounts[my_rank], sendtypes[my_rank],
                               my_rank, tag, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);

        if ((mpierr = MPI_Wait(rcvids, &status)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
#endif
    }

    /* When send to self is complete there is nothing left to do if
     * ntasks==1. */
    if (ntasks == 1)
        return PIO_NOERR;

    for (int i = 0; i < ntasks; i++)
    {
        rcvids[i] = MPI_REQUEST_NULL;
        swapids[i] = 0;
    }

    if (isend)
        for (int i = 0; i < ntasks; i++)
            sndids[i] = MPI_REQUEST_NULL;

    if (handshake)
        for (int i = 0; i < ntasks; i++)
            hs_rcvids[i] = MPI_REQUEST_NULL;

    steps = 0;
    for (istep = 0; istep < ceil2(ntasks) - 1; istep++)
    {
        p = pair(ntasks, istep, my_rank);
        if (p >= 0 && (sendcounts[p] > 0 || recvcounts[p] > 0))
            swapids[steps++] = p;
    }

    if (steps == 1)
    {
        maxreq = 1;
        maxreqh = 1;
    }
    else
    {
        if (max_requests > 1 && max_requests < steps)
        {
            maxreq = max_requests;
            maxreqh = maxreq / 2;
        }
        else if (max_requests >= steps)
        {
            maxreq = steps;
            maxreqh = steps;
        }
        else
        {
            maxreq = 2;
            maxreqh = 1;
        }
    }

    /* If handshaking is in use, do a nonblocking recieve to listen
     * for it. */
    if (handshake)
    {
        for (istep = 0; istep < maxreq; istep++)
        {
            p = swapids[istep];
            if (sendcounts[p] > 0)
            {
                tag = my_rank + offset_t;
                if ((mpierr = MPI_Irecv(&hs, 1, MPI_INT, p, tag, comm, hs_rcvids + istep)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            }
        }
    }

    /* Post up to maxreq irecv's. */
    for (istep = 0; istep < maxreq; istep++)
    {
        p = swapids[istep];
        if (recvcounts[p] > 0)
        {
            tag = p + offset_t;
            ptr = (char *)recvbuf + rdispls[p];

            if ((mpierr = MPI_Irecv(ptr, recvcounts[p], recvtypes[p], p, tag, comm,
                                    rcvids + istep)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);

            if (handshake)
                if ((mpierr = MPI_Send(&hs, 1, MPI_INT, p, tag, comm)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        }
    }

    /* Tell the paired task that this tasks' has posted it's irecvs'. */
    rstep = maxreq;
    for (istep = 0; istep < steps; istep++)
    {
        p = swapids[istep];
        if (sendcounts[p] > 0)
        {
            tag = my_rank + offset_t;
            /* If handshake is enabled don't post sends until the
             * receiving task has posted recvs. */
            if (handshake)
            {
                if ((mpierr = MPI_Wait(hs_rcvids + istep, &status)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                hs_rcvids[istep] = MPI_REQUEST_NULL;
            }
            ptr = (char *)sendbuf + sdispls[p];

            if (isend)
            {                
                if ((mpierr = MPI_Irsend(ptr, sendcounts[p], sendtypes[p], p, tag, comm,
                                         sndids + istep)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            }
            else
            {
                if ((mpierr = MPI_Send(ptr, sendcounts[p], sendtypes[p], p, tag, comm)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            }
        }

        /* We did comms in sets of size max_reqs, if istep > maxreqh
         * then there is a remainder that must be handled. */
        if (istep > maxreqh)
        {
            p = istep - maxreqh;
            if (rcvids[p] != MPI_REQUEST_NULL)
            {
                if ((mpierr = MPI_Wait(rcvids + p, &status)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                rcvids[p] = MPI_REQUEST_NULL;
            }
            if (rstep < steps)
            {
                p = swapids[rstep];
                if (handshake && sendcounts[p] > 0)
                {
                    tag = my_rank + offset_t;
                    if ((mpierr = MPI_Irecv(&hs, 1, MPI_INT, p, tag, comm, hs_rcvids+rstep)))
                        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                }
                if (recvcounts[p] > 0)
                {
                    tag = p + offset_t;

                    ptr = (char *)recvbuf + rdispls[p];
                    if ((mpierr = MPI_Irecv(ptr, recvcounts[p], recvtypes[p], p, tag, comm, rcvids + rstep)))
                        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                    if (handshake)
                        if ((mpierr = MPI_Send(&hs, 1, MPI_INT, p, tag, comm)))
                            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                }
                rstep++;
            }
        }
    }

    /* If steps > 0 there are still outstanding messages, wait for
     * them here. */
    if (steps > 0)
    {
        if ((mpierr = MPI_Waitall(steps, rcvids, MPI_STATUSES_IGNORE)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if (isend)
            if ((mpierr = MPI_Waitall(steps, sndids, MPI_STATUSES_IGNORE)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    }

    return PIO_NOERR;
}

/**
 * Provides the functionality of MPI_Gatherv with flow control
 * options. This function is not currently used, but we hope it will
 * be useful in future optimizations.
 *
 * @param sendbuf starting address of send buffer.
 * @param sendcnt number of elements in send buffer.
 * @param sendtype data type of send buffer elements.
 * @param recvbuf address of receive buffer.
 * @param recvcnts integer array (of length group size) containing the
 * number of elements that are received from each process (significant
 * only at root).
 * @param displs integer array (of length group size). Entry i
 * specifies the displacement relative to recvbuf at which to place
 * the incoming data from process i (significant only at root).
 * @param recvtype data type of recv buffer elements (significant only
 * at root).
 * @param root rank of receiving process.
 * @param comm communicator.
 * @param flow_cntl if non-zero, flow control will be used.
 * @returns 0 for success, error code otherwise.
 */
int pio_fc_gatherv(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                   void *recvbuf, const int *recvcnts, const int *displs,
                   MPI_Datatype recvtype, int root, MPI_Comm comm, int flow_cntl)
{
    bool fc_gather;
    int gather_block_size;
    int mytask, nprocs;
    int mtag;
    MPI_Status status;
    int hs;
    int dsize;
    int mpierr;  /* Return code from MPI functions. */

    if (flow_cntl > 0)
    {
        fc_gather = true;
        gather_block_size = min(flow_cntl, MAX_GATHER_BLOCK_SIZE);
    }
    else
    {
        fc_gather = false;
    }

    if (fc_gather)
    {
        if ((mpierr = MPI_Comm_rank(comm, &mytask)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        if ((mpierr = MPI_Comm_size(comm, &nprocs)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);

        mtag = 2 * nprocs;
        hs = 1;

        if (mytask == root)
        {
            int preposts = min(nprocs-1, gather_block_size);
            int head = 0;
            int count = 0;
            int tail = 0;
            MPI_Request rcvid[gather_block_size];

            if ((mpierr = MPI_Type_size(recvtype, &dsize)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);

            for (int p = 0; p < nprocs; p++)
            {
                if (p != root)
                {
                    if (recvcnts[p] > 0)
                    {
                        count++;
                        if (count > preposts)
                        {
                            if ((mpierr = MPI_Wait(rcvid + tail, &status)))
                                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                            tail = (tail + 1) % preposts;
                        }

                        void *ptr = (void *)((char *)recvbuf + dsize * displs[p]);

                        if ((mpierr = MPI_Irecv(ptr, recvcnts[p], recvtype, p, mtag, comm, rcvid + head)))
                            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                        head = (head + 1) % preposts;
                        if ((mpierr = MPI_Send(&hs, 1, MPI_INT, p, mtag, comm)))
                            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                    }
                }
            }

            /* copy local data */
            if ((mpierr = MPI_Type_size(sendtype, &dsize)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            if ((mpierr = MPI_Sendrecv(sendbuf, sendcnt, sendtype, mytask, 102, recvbuf, recvcnts[mytask],
                                       recvtype, mytask, 102, comm, &status)))
                return check_mpi(NULL, mpierr, __FILE__, __LINE__);

            count = min(count, preposts);
            if (count > 0)
                if ((mpierr = MPI_Waitall(count, rcvid, MPI_STATUSES_IGNORE)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
        }
        else
        {
            if (sendcnt > 0)
            {
                if ((mpierr = MPI_Recv(&hs, 1, MPI_INT, root, mtag, comm, &status)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
                if ((mpierr = MPI_Send(sendbuf, sendcnt, sendtype, root, mtag, comm)))
                    return check_mpi(NULL, mpierr, __FILE__, __LINE__);
            }
        }
    }
    else
    {
        if ((mpierr = MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts,
                                  displs, recvtype, root, comm)))
            return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    }

    return PIO_NOERR;
}

