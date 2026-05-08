struct RK2Workspace {
    std::vector<double> kX1, kX2, qX;
    void ensure_size(std::size_t n) {
        kX1.resize(n);
        kX2.resize(n);
        qX.resize(n);
    }
};

void RK2( vector<double>& X, double& pt, RK2Workspace& rk2wk)
{
    rk2wk.ensure_size(X.size());
    const std::size_t n = X.size();
    onestep(X, rk2wk.kX1, pt );   
    for(size_t i = 0; i < n; ++i) rk2wk.qX[i] = X[i] + Dt * rk2wk.kX1[i];    
    onestep(rk2wk.qX, rk2wk.kX2, pt + HDt );
    for(size_t i = 0; i < n; ++i) X[i] += HDt * (rk2wk.kX1[i] + rk2wk.kX2[i]);
    pt += Dt;
}
