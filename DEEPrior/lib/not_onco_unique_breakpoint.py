def not_onco_unique_breakpoint(not_onco):
    """
    :param not_onco: huge NotOnco list
    :return:pruned Not onco list with only unique breakpoint. If multiple fusions have the same breakpoint, tissue is
    updated with the names of all the tissue involved
    """
    not_onco_pruned_obj = []
    fusion_list = []

    for i in not_onco:
        if i.tissue not in ['ESC', 'MSC', 'MFC10']:
            fus = "%s_%s_%s_%s" % (i.chr5p, i.coord5p, i.chr3p, i.coord3p)
            if fus not in fusion_list:
                fusion_list.append(fus)
                not_onco_pruned_obj.append(i)
            else:
                ind = fusion_list.index(fus)
                not_onco_pruned_obj[ind].tissue += '_%s' % i.tissue

    return not_onco_pruned_obj
